#include <ros/ros.h>
#include <ros/package.h>
#include <ros/console.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <chrono>
#include <random>
#include <string>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <sensor_msgs/PointCloud2.h>

#include <decomp_ros_utils/data_ros_utils.h>
#include <decomp_util/ellipsoid_decomp.h>
#include "mapping/geoutils.hpp"

#include "trajectory_utils/poly_traj_utils.hpp"
#include "traj_opt/global_traj_opt.h"

#include "wr_msg/wr_msg.hpp"
#include "visualization/visualization.hpp"
#include <std_msgs/Float32.h>

#include <quadrotor_msgs/PolynomialTrajectory.h>
#include <quadrotor_msgs/HpolyTimeSerial.h>
#include "corridor_helper/corridor_helper.h"

enum PCD_TYPE
{
  PCD = 0,
  MSG = 1
};

pcl::KdTreeFLANN<pcl::PointXYZ> kdtreeLocalMap;
std::vector<int> pointIdxRadiusSearch;
std::vector<float> pointRadiusSquaredDistance;

// generated
std::shared_ptr<traj_opt::GlobTrajOpt> globTrajOptPtr_;
traj_opt::Config config_;
Trajectory global_traj_;
ros::Time glb_start_timestamp_;

Eigen::Vector3d target_p_, initial_p_;
std::vector<Eigen::Vector3d> gate_pos_list_;

std::vector<Eigen::MatrixXd> hPolys_;
std::vector<Eigen::Vector3d> path_;

double inflate_size_ = 1;

std::shared_ptr<visualization::Visualization> visPtr_;

double l_corridor_;
double bbox_width_;
double bbox_width_x_, bbox_width_y_, bbox_width_z_;
double desired_len_;

ros::Publisher all_map_pub_, hPolyPub_, gridmap_vis_pub_, gridmap_inflate_vis_pub_;
ros::Subscriber looking_sub;

double resolution_;
double x_size, y_size, z_size;
double z_floor_, z_ceil_;
double x_ceil_, x_floor_;
double y_ceil_, y_floor_;

double drone_radius_;

pcl::PointCloud<pcl::PointXYZ> cloudDense;
std::shared_ptr<corridor_opt::CorridorHelper> CorridorPtr_;
std::vector<Eigen::Matrix3d> best_time_serial_;
std::vector<int> gate_signs;

int looking_id_ = 0;
quadrotor_msgs::HpolyTimeSerial HpolyTimeserial;
double extend_d_backward_, extend_d_forward_;

void compressPoly(Polyhedron3D &poly, double dx)
{
  vec_E<Hyperplane3D> hyper_planes = poly.hyperplanes();
  for (uint j = 0; j < hyper_planes.size(); j++)
  {
    hyper_planes[j].p_ = hyper_planes[j].p_ - hyper_planes[j].n_ * dx;
  }
  poly = Polyhedron3D(hyper_planes);
}

void compressPoly(Eigen::MatrixXd &poly, double dx)
{
  for (int i = 0; i < poly.cols(); ++i)
  {
    poly.col(i).tail(3) = poly.col(i).tail(3) - poly.col(i).head(3) * dx;
  }
}

void getPointCloudAroundLineWithKDtree(const vec_Vec3f &line, const double maxWidth, vec_Vec3f &pc)
{
  pc.clear();
  Eigen::Vector3d p0 = line.front();
  Eigen::Vector3d p1 = line.back();
  Eigen::Vector3d p_center = 0.5 * (p0 + p1);
  double sensing_range = (p_center - p0).norm() + maxWidth;
  Eigen::Vector3d tmp_p;

  std::vector<int> pointIdxRadiusSearch;
  std::vector<float> pointRadiusSquaredDistance;

  pcl::PointXYZ searchPoint(p_center[0], p_center[1], p_center[2]);
  pcl::PointXYZ pt;
  if (kdtreeLocalMap.radiusSearch(searchPoint, sensing_range,
                                  pointIdxRadiusSearch,
                                  pointRadiusSquaredDistance) > 0)
  {
    for (size_t i = 0; i < pointIdxRadiusSearch.size(); ++i)
    {
      tmp_p << cloudDense.points[pointIdxRadiusSearch[i]].x, cloudDense.points[pointIdxRadiusSearch[i]].y, cloudDense.points[pointIdxRadiusSearch[i]].z;
      pc.push_back(tmp_p);
    }
  }

  double search_step = 0.5;
  for (tmp_p.x() = p_center.x() - sensing_range; tmp_p.x() <= p_center.x() + sensing_range; tmp_p.x() += search_step)
  {
    for (tmp_p.y() = p_center.y() - sensing_range; tmp_p.y() <= p_center.y() + sensing_range; tmp_p.y() += search_step)
    {
      for (tmp_p.z() = p_center.z() - sensing_range; tmp_p.z() <= p_center.z() + sensing_range; tmp_p.z() += search_step)
      {
        if (tmp_p.x() <= x_floor_ || tmp_p.x() >= x_ceil_)
        {
          pc.push_back(tmp_p);
          continue;
        }
        if (tmp_p.y() <= y_floor_ || tmp_p.y() >= y_ceil_)
        {
          pc.push_back(tmp_p);
          continue;
        }
        if (tmp_p.z() <= z_floor_ || tmp_p.z() >= z_ceil_)
        {
          pc.push_back(tmp_p);
          continue;
        }
      }
    }
  }
}

void SFCgenerate(std::vector<Eigen::MatrixXd> &hPolys, std::vector<Eigen::Vector3d> &path)
{
  double real_len;
  Eigen::Vector3d pos_backward, pos_forward, forward_dir;
  Eigen::VectorXd duration = global_traj_.getDurations();
  int num = duration.size();
  double accumulate_t = 0.0, accumulate_dist = 0.0;
  double dt = 0.01;
  std::vector<std::pair<double, double>> time_pair_set;

  // int current_id = 0;
  for (int i = 0; i < num; i++)
  {
    double dist = global_traj_.getDistance(0.01, accumulate_t, accumulate_t + duration[i]);
    int N = dist / desired_len_;
    N = N > 0 ? N : 1;
    std::cout << "gate_index: " << i << " N: " << N << std::endl;
    real_len = dist / N;
    std::cout << "dist: " << dist << " real_len: " << real_len << std::endl;

    Eigen::Vector3d pos_start = global_traj_.getPos(accumulate_t);

    if (i == 0)
    {
      path.push_back(pos_start);
      gate_signs.push_back(0);
    }
    else
    {
      path.push_back(pos_start - 0.2 * global_traj_.getVel(accumulate_t).normalized());
      gate_signs.push_back(0);
      path.push_back(pos_start);
      gate_signs.push_back(1);

      path.push_back(pos_start + 0.2 * global_traj_.getVel(accumulate_t).normalized());
      gate_signs.push_back(0);
    }

    Eigen::Vector3d p;
    accumulate_dist = 0.0;
    double dist_tmp = 0.0;
    int n = 0;

    for (double t = dt; t < duration[i]; t += dt)
    {
      if (n == N - 1)
        break;
      p = global_traj_.getPos(accumulate_t + t);
      accumulate_dist += (p - pos_start).norm();
      if ((accumulate_dist - dist_tmp) >= real_len)
      {
        path.push_back(p);
        gate_signs.push_back(0);
        dist_tmp += real_len;
        n++;
      }
      pos_start = p;
    }
    accumulate_t += duration[i];
  }

  path.push_back(target_p_);

  for (int i = 0; i < 100; i++)
  {
    visPtr_->visualize_pointcloud(path, "sampled_path_points");
    ros::Duration(0.01).sleep();
  }

  Eigen::Vector3d p_last = path[0];
  for (int i = 1; i < int(path.size()); i++)
  {
    p_last = path[i];
  }

  vec_Vec3f obs_pc;
  EllipsoidDecomp3D decomp_util;
  std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> keyPts;
  decomp_util.set_local_bbox(Eigen::Vector3d(bbox_width_x_, bbox_width_y_, bbox_width_z_));

  double maxWidth = std::max(std::max(bbox_width_x_, bbox_width_y_), bbox_width_z_);

  vec_E<Polyhedron3D> decompPolys;

  int path_len = path.size();

  int idx = 0;
  keyPts.clear();
  while (idx < path_len - 1)
  {
    int next_idx = idx + 1;

    // generate corridor with idx and next_idx
    vec_Vec3f line;
    line.push_back(path[idx]);
    line.push_back(path[next_idx]);
    keyPts.emplace_back(path[idx], path[next_idx]);
    getPointCloudAroundLineWithKDtree(line, maxWidth, obs_pc);
    decomp_util.set_obs(obs_pc);
    decomp_util.dilate(line);
    Polyhedron3D poly = decomp_util.get_polyhedrons()[0];
    decompPolys.push_back(poly);

    // find a farest idx in current corridor
    // idx ++;
    idx = next_idx;
    //   while (idx + 1 < path_len && decompPolys.back().inside(path[idx + 1])) {
    //     idx++;
    //   }
  }

  hPolys.clear();
  Eigen::MatrixXd current_poly;
  for (uint i = 0; i < decompPolys.size(); i++)
  {
    vec_E<Hyperplane3D> current_hyperplanes = decompPolys[i].hyperplanes();
    current_poly.resize(6, current_hyperplanes.size());
    for (uint j = 0; j < current_hyperplanes.size(); j++)
    {
      current_poly.col(j) << current_hyperplanes[j].n_, current_hyperplanes[j].p_;
      // outside
    }
    hPolys.push_back(current_poly);
  }

  //------------- directly compress the corridor with the radius of the quadrotor --------------------
  Eigen::MatrixXd curIH;
  Eigen::Vector3d interior;
  for (int i = 0; i < (int)hPolys.size(); i++)
  {
    std::cout << i << "-th hpoly: " << hPolys[i].cols() << std::endl;
    if (geoutils::findInteriorDist(current_poly, interior) < drone_radius_)
    {
      std::cout << "feasible interior dist small than " << drone_radius_ << " m" << std::endl;
    }
    else
    {
      compressPoly(hPolys[i], drone_radius_);
    }
  }
}

void global_plan()
{
  Eigen::Vector3d start_p, start_v, start_a;
  Eigen::Vector3d target_p, target_v, target_a;
  Eigen::MatrixXd initState(3, 3), finState(3, 3);

  start_p = initial_p_;
  start_v << 0, 0, 0;
  start_a << 0, 0, 0;

  target_p = target_p_;
  target_v << 0, 0, 0;
  target_a << 0, 0, 0;

  initState << start_p, start_v, start_a;
  finState << target_p, target_v, target_a;

  // N is piece num of global traj
  int N = gate_pos_list_.size() + 1;

  // waypoints include start_p and target_p
  Eigen::MatrixXd waypoints;
  Eigen::VectorXd Ts;

  waypoints.resize(3, N + 1);
  waypoints.col(0) = start_p;
  for (int i = 0; i < N - 1; i++)
  {
    waypoints.col(i + 1) = gate_pos_list_[i];
  }

  waypoints.col(N) = target_p;

  Ts.resize(N);
  for (int i = 0; i < N; ++i)
  {
    Ts[i] = (waypoints.col(i + 1) - waypoints.col(i)).norm() / config_.max_vel;
  }

  std::cout << "start global plan !" << std::endl;
  bool res = globTrajOptPtr_->generate_traj(initState, finState, waypoints.block(0, 1, 3, N - 1), Ts, global_traj_);

  if (res)
  {
    glb_start_timestamp_ = ros::Time::now();
  }
}

void corridor_assignment(sensor_msgs::PointCloud2 &time_point)
{
  pcl::PointXYZ pt;
  pcl::PointCloud<pcl::PointXYZ> pcd;
  int idx_last = -1;
  double accumulate_dist = 0;
  int cnt = 0;
  Eigen::Vector3d pos_last = global_traj_.getPos(0);

  for (double t = 0; t <= global_traj_.getTotalDuration(); t += 0.01)
  {
    Eigen::Vector3d pos = global_traj_.getPos(t);
    int id_tmp = idx_last < 0 ? 0 : idx_last;
    int idx_best = CorridorPtr_->getBestHpoly(pos, id_tmp, 10);
    if (idx_best > idx_last)
    {
      pt.x = pos.x();
      pt.y = pos.y();
      pt.z = pos.z();
      pcd.push_back(pt);

      best_time_serial_[cnt].col(0)[0] = t;
      best_time_serial_[cnt].col(1) = pos;
      if (cnt > 0)
      {
        best_time_serial_[cnt - 1].col(0)[1] = t;
        best_time_serial_[cnt - 1].col(0)[2] = accumulate_dist;
        best_time_serial_[cnt - 1].col(2) = pos;
        accumulate_dist = 0;
      }
      cnt++;
      idx_last++;
    }
    accumulate_dist += (pos - pos_last).norm();
    pos_last = pos;
  }

  best_time_serial_[cnt - 1].col(0)[1] = global_traj_.getTotalDuration();
  best_time_serial_[cnt - 1].col(0)[2] = accumulate_dist;
  best_time_serial_[cnt - 1].col(2) = global_traj_.getPos(global_traj_.getTotalDuration());

  pcl::toROSMsg(pcd, time_point);
  time_point.header.frame_id = "world";
}

quadrotor_msgs::PolynomialTrajectory traj2msg(const Trajectory &traj, const ros::Time &traj_timestamp)
{
  static int count = 0;
  quadrotor_msgs::PolynomialTrajectory traj_msg;
  traj_msg.header.seq = count;
  traj_msg.header.stamp = traj_timestamp;
  traj_msg.trajectory_id = count;
  traj_msg.action = quadrotor_msgs::PolynomialTrajectory::ACTION_ADD;
  traj_msg.num_order = traj[0].getOrder(); // the order of polynomial
  traj_msg.num_segment = traj.getPieceNum();
  traj_msg.start_yaw = 0;
  traj_msg.final_yaw = 0;
  for (unsigned int i = 0; i < traj_msg.num_segment; i++)
  {
    for (int j = 0; j <= traj[i].getOrder(); j++)
    {
      //   CoefficientMat coemat = traj[i].getCoeffMat(true);
      CoefficientMat coemat = traj[i].getCoeffMat(false);
      traj_msg.coef_x.push_back(coemat(0, j));
      traj_msg.coef_y.push_back(coemat(1, j));
      traj_msg.coef_z.push_back(coemat(2, j));
    }
    traj_msg.time.push_back(traj[i].getDuration());
    traj_msg.order.push_back(traj[i].getOrder());
  }
  traj_msg.mag_coeff = 1;
  count++;
  return traj_msg;
}

inline decomp_ros_msgs::PolyhedronArray HpolysToMsg(const std::vector<Eigen::MatrixXd> &hPolys, ros::Time time_stamp)
{
  vec_E<Polyhedron3D> decompPolys;
  decomp_ros_msgs::PolyhedronArray poly_msg;
  for (const auto &poly : hPolys)
  {
    vec_E<Hyperplane3D> hyper_planes;
    hyper_planes.resize(poly.cols());
    for (int i = 0; i < poly.cols(); ++i)
    {
      hyper_planes[i].n_ = poly.col(i).head(3);
      hyper_planes[i].p_ = poly.col(i).tail(3);
    }
    decompPolys.emplace_back(hyper_planes);
  }

  poly_msg = DecompROS::polyhedron_array_to_ros(decompPolys);
  poly_msg.header.frame_id = "world";
  poly_msg.header.stamp = time_stamp;

  return poly_msg;
}

inline void visCorridor(const vec_E<Polyhedron3D> &polyhedra)
{
  decomp_ros_msgs::PolyhedronArray poly_msg = DecompROS::polyhedron_array_to_ros(polyhedra);
  poly_msg.header.frame_id = "world";
  poly_msg.header.stamp = ros::Time::now();
  hPolyPub_.publish(poly_msg);
}

inline void visCorridor(const std::vector<Eigen::MatrixXd> &hPolys)
{
  vec_E<Polyhedron3D> decompPolys;
  for (const auto &poly : hPolys)
  {
    vec_E<Hyperplane3D> hyper_planes;
    hyper_planes.resize(poly.cols());
    for (int i = 0; i < poly.cols(); ++i)
    {
      hyper_planes[i].n_ = poly.col(i).head(3);
      hyper_planes[i].p_ = poly.col(i).tail(3);
    }
    decompPolys.emplace_back(hyper_planes);
  }
  visCorridor(decompPolys);
}

void lookingIdCallback(const std_msgs::Float32 msg)
{
  looking_id_ = msg.data;
  ROS_WARN("looking_id_ update !");
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "random_map_sensing");
  ros::NodeHandle n("~");

  looking_sub = n.subscribe("looking_data", 1, lookingIdCallback);

  all_map_pub_ = n.advertise<sensor_msgs::PointCloud2>("global_cloud", 1);
  hPolyPub_ = n.advertise<decomp_ros_msgs::PolyhedronArray>("polyhedra", 1);
  gridmap_vis_pub_ = n.advertise<sensor_msgs::PointCloud2>("vis_gridmap", 1);
  gridmap_inflate_vis_pub_ = n.advertise<sensor_msgs::PointCloud2>("inflate_vis_gridmap", 1);

  CorridorPtr_ = std::make_shared<corridor_opt::CorridorHelper>(n);

  double inflate_size;
  n.param("map/x_size", x_size, 50.0);
  n.param("map/y_size", y_size, 50.0);
  n.param("map/z_size", z_size, 5.0);
  n.param("map/resolution", resolution_, 0.1);
  n.param("map/inflate_size", inflate_size, 0.3);

  std::vector<double> target_p_vec;
  n.getParam("target_pos", target_p_vec);
  target_p_ << target_p_vec[0], target_p_vec[1], target_p_vec[2];

  std::vector<double> init_p_vec;
  n.getParam("initial_pos", init_p_vec);
  initial_p_ << init_p_vec[0], init_p_vec[1], init_p_vec[2];

  std::vector<double> tmp_gate_pos_vec;
  n.getParam("gate_pos_vec", tmp_gate_pos_vec);
  int num = tmp_gate_pos_vec.size() / 3;
  gate_pos_list_.resize(num);
  for (int i = 0; i < num; i++)
  {
    gate_pos_list_[i] << tmp_gate_pos_vec[3 * i], tmp_gate_pos_vec[3 * i + 1], tmp_gate_pos_vec[3 * i + 2];
  }

  n.param("z_ceil", z_ceil_, 1.0);
  n.param("z_floor", z_floor_, 0.0);
  n.param("x_ceil", x_ceil_, 2.0);
  n.param("x_floor", x_floor_, -2.0);
  n.param("y_ceil", y_ceil_, 2.0);
  n.param("y_floor", y_floor_, -2.0);

  n.param("corridor/l_corridor", l_corridor_, 2.0);
  n.param("corridor/bbox_width_x", bbox_width_x_, 4.0);
  n.param("corridor/bbox_width_y", bbox_width_y_, 4.0);
  n.param("corridor/bbox_width_z", bbox_width_z_, 4.0);

  n.param("corridor/desired_len", desired_len_, 2.0);
  n.param("corridor/drone_radius", drone_radius_, 0.5);
  n.param("corridor/extend_d_backward", extend_d_backward_, 0.2);
  n.param("corridor/extend_d_forward", extend_d_forward_, 0.2);

  visPtr_ = std::make_shared<visualization::Visualization>(n);

  config_.load(n);
  globTrajOptPtr_ = std::make_shared<traj_opt::GlobTrajOpt>(n, config_);
  globTrajOptPtr_->setVisualizer(visPtr_);

  ROS_WARN("start loading pointclouds");
  sensor_msgs::PointCloud2 globalMap_pcd;

  int loader_type = PCD_TYPE::MSG;
  n.param("pcd_type", loader_type, 1);
  if (loader_type == PCD_TYPE::PCD)
  {
    //--------------------------------------------------- pcd loader -----------------------------------------------------------------
    std::string path_to_pcd = ros::package::getPath("plan_manage") + "/pcd/globalmap.pcd";
    pcl::io::loadPCDFile<pcl::PointXYZ>(path_to_pcd, cloudDense);

    cloudDense.width = cloudDense.size();
    cloudDense.height = 1;
    cloudDense.is_dense = true;
  }
  else if (loader_type == PCD_TYPE::MSG)
  {
    //-------------------------------------------------- msg loader ------------------------------------------------------------------
    wr_msg::readMsg(globalMap_pcd, ros::package::getPath("plan_manage") + "/offlines/globalmap.bin");
    pcl::fromROSMsg(globalMap_pcd, cloudDense);
  }
  else
  {
    ROS_ERROR("Invalid pcd loader type !");
    return 0;
  }

  ROS_WARN("finish loading pointclouds");
  ROS_WARN("start setting map");

  kdtreeLocalMap.setInputCloud(cloudDense.makeShared());

  ROS_WARN("finish setting map");

  ROS_WARN("Start generating global trajectory");
  global_plan();
  ROS_WARN("Finish generating global trajectory");

  std::cout << "max_vel: " << global_traj_.getMaxVelRate() << std::endl;

  ROS_WARN("Start generating SFC");
  SFCgenerate(hPolys_, path_);
  ROS_WARN("Finish generating SFC");

  CorridorPtr_->setCorridor(hPolys_);
  best_time_serial_.resize(hPolys_.size());
  std::fill(best_time_serial_.begin(), best_time_serial_.end(), Eigen::Matrix3d::Zero());

  sensor_msgs::PointCloud2 time_point;
  corridor_assignment(time_point);

  for (int i = 0; i < int(hPolys_.size()); i++)
  {
    if (gate_signs[i] == 1)
    {
      best_time_serial_[i].col(0)[0] = best_time_serial_[i].col(0)[0] - extend_d_backward_;
      best_time_serial_[i].col(1) = global_traj_.getPos(best_time_serial_[i].col(0)[0]);
      for (int j = i; j > 0; j--)
      {
        if (best_time_serial_[j - 1].col(0)[1] > best_time_serial_[i].col(0)[0])
        {
          best_time_serial_[j - 1].col(0)[1] = best_time_serial_[i].col(0)[0];
          best_time_serial_[j - 1].col(2) = global_traj_.getPos(best_time_serial_[i].col(0)[0]);
        }
        else
          break;
      }
    }
  }

  for (int i = 0; i < int(hPolys_.size()); i++)
  {
    if (gate_signs[i] == 1)
    {
      best_time_serial_[i].col(0)[1] = best_time_serial_[i].col(0)[1] + extend_d_forward_;
      best_time_serial_[i].col(2) = global_traj_.getPos(best_time_serial_[i].col(0)[1]);
      for (int j = i; j < int(hPolys_.size()) - 1; j++)
      {
        if (best_time_serial_[j + 1].col(0)[0] < best_time_serial_[i].col(0)[1])
        {
          best_time_serial_[j + 1].col(0)[0] = best_time_serial_[i].col(0)[1];
          best_time_serial_[j + 1].col(1) = global_traj_.getPos(best_time_serial_[i].col(0)[1]);
        }
        else
          break;
      }
    }
  }

  std::vector<Eigen::MatrixXd> hpolys_tmp;
  std::vector<Eigen::Matrix3d> best_time_serial_tmp_;
  int numh = hPolys_.size();
  for (int i = 0; i < numh; i++)
  {
    if (best_time_serial_[i].col(0)[1] < best_time_serial_[i].col(0)[0])
      continue;
    hpolys_tmp.push_back(hPolys_[i]);
    best_time_serial_tmp_.push_back(best_time_serial_[i]);
  }

  hPolys_ = hpolys_tmp;
  best_time_serial_ = best_time_serial_tmp_;

  double dist_rr = 0;
  std::cout << "hpolys time searial: " << std::endl;
  for (int i = 0; i < int(hPolys_.size()); i++)
  {
    std::cout << "hpoly " << i << " : t_start: " << best_time_serial_[i].col(0)[0] << " t_end: " << best_time_serial_[i].col(0)[1] << " dist: " << best_time_serial_[i].col(0)[2] << std::endl;
    dist_rr += best_time_serial_[i].col(0)[2];
    HpolyTimeserial.t_start.push_back(best_time_serial_[i].col(0)[0]);
    HpolyTimeserial.t_end.push_back(best_time_serial_[i].col(0)[1]);
  }

  std::cout << "total dist calculated: " << dist_rr << std::endl;
  std::cout << "total dist ground truth: " << global_traj_.getTotalDistance(0.01) << std::endl;

  ROS_WARN("Start writing msgs to binary");
  // write the traj and corridor msg to the binary
  quadrotor_msgs::PolynomialTrajectory traj_msg_for_write;
  traj_msg_for_write = traj2msg(global_traj_, glb_start_timestamp_);
  wr_msg::writeMsg(traj_msg_for_write, ros::package::getPath("plan_manage") + "/offlines/global_trajectory.bin");

  decomp_ros_msgs::PolyhedronArray poly_msg_for_write;
  poly_msg_for_write = HpolysToMsg(hPolys_, ros::Time::now());
  wr_msg::writeMsg(poly_msg_for_write, ros::package::getPath("plan_manage") + "/offlines/corridors.bin");

  wr_msg::writeMsg(HpolyTimeserial, ros::package::getPath("plan_manage") + "/offlines/HpolyTimeSerial.bin");
  ROS_WARN("Finish writing msgs to binary");

  ros::Rate loop_rate(20);

  ROS_WARN("visualization start!");

  ROS_WARN("Start writing map clouds to binary");
  pcl::toROSMsg(cloudDense, globalMap_pcd);
  globalMap_pcd.header.frame_id = "world";
  wr_msg::writeMsg(globalMap_pcd, ros::package::getPath("plan_manage") + "/offlines/globalmap.bin");
  ROS_WARN("Finish writing map clouds to binary");

  // visualization for test.
  while (ros::ok())
  {

    visPtr_->visualize_traj(global_traj_, "global_traj_generated");

    static ros::Time t1 = ros::Time::now();
    if ((ros::Time::now() - t1).toSec() > 2.0)
    {
      t1 = ros::Time::now();
      gridmap_vis_pub_.publish(globalMap_pcd);
      all_map_pub_.publish(time_point);
      visCorridor(hPolys_);
    }

    visPtr_->visualize_path(path_, "path_for_polyhedra");
    visPtr_->visualize_pointcloud(path_, "sampled_path_points");

    ros::spinOnce();
    loop_rate.sleep();
  }
}