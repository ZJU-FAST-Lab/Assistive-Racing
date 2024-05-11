#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
// #include <pcl/search/kdtree.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl_conversions/pcl_conversions.h>
#include <iostream>

#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/Vector3.h>
#include <math.h>
#include <nav_msgs/Odometry.h>
#include <ros/console.h>
#include <ros/ros.h>
#include <ros/package.h>
#include <sensor_msgs/PointCloud2.h>
#include <quadrotor_msgs/PolynomialTrajectory.h>
#include <Eigen/Eigen>
#include <random>

#include <decomp_ros_utils/data_ros_utils.h>
#include <decomp_util/ellipsoid_decomp.h>
#include "visualization/visualization.hpp"

#include "trajectory_utils/poly_traj_utils.hpp"
#include "traj_opt/global_traj_opt.h"

#include "wr_msg/wr_msg.hpp"

#include <octomap/octomap.h>

using namespace std;

// pcl::search::KdTree<pcl::PointXYZ> kdtreeLocalMap;
pcl::KdTreeFLANN<pcl::PointXYZ> kdtreeLocalMap;
vector<int> pointIdxRadiusSearch;
vector<float> pointRadiusSquaredDistance;

random_device rd;
default_random_engine eng(rd());
uniform_real_distribution<double> rand_x;
uniform_real_distribution<double> rand_y;
uniform_real_distribution<double> rand_w;
uniform_real_distribution<double> rand_h;

ros::Publisher _all_map_pub;
ros::Publisher gridmap_vs_pub, gridmap_inflate_vs_pub;
ros::Publisher click_map_pub_;
ros::Publisher global_traj_pub_;
ros::Subscriber _odom_sub;
// ros::Subscriber global_traj_sub_;
ros::Publisher circle_pub_;
ros::Time t_tic;
ros::Publisher hPolyPub_;

vector<double> _state;

int _obs_num;
double _x_size, _y_size, _z_size;
double _sensing_range, _resolution, _sense_rate;
double _init_x, _init_y, _init_z;
double inner_gate_radius_, outter_gate_radius_;

bool _map_ok = false;
bool _has_odom = false;

int circle_num_;
double radius_l_, radius_h_, z_l_, z_h_;
double theta_;

// Gates:[x, y, z, theta, signal]
Eigen::MatrixXd Gates;
std::vector<double> gate_vectors_;
int idx_now = 0;

sensor_msgs::PointCloud2 globalMap_pcd;
pcl::PointCloud<pcl::PointXYZ> cloudMap;
pcl::PointCloud<pcl::PointXYZ> clicked_cloud_;

// Grid map related
double inflate_size_ = 1;

// corridor related
double l_corridor_ = 1;
double bbox_width_ = 3;
std::shared_ptr<visualization::Visualization> visPtr_;

Trajectory global_traj_;
bool if_SFC_generated_ = false;
bool get_glb_traj_ = false;

std::vector<Eigen::MatrixXd> hPolys_;
std::vector<Eigen::Vector3d> path_;

std::vector<double> target_p_vec_;
Eigen::Vector3d target_p_;
ros::Time glb_start_timestamp_;

double desired_len_;
bool remove_floor_ceil_ = false;

std::vector<bool> fake_gates_;

// Title

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

void column_genreate(const Eigen::Vector3d pos, const double inflate_r)
{
  pcl::PointXYZ racingPt;
  Eigen::Vector3d Cpt;
  double inflate_num = inflate_r / _resolution;
  Cpt(0) = pos.x();
  Cpt(1) = pos.y();
  Cpt(2) = 0.0;
  for (double zm = 0; zm < pos.z(); zm += _resolution)
  {
    Cpt(2) = zm;
    Eigen::Vector3d Cpt_if;
    for (int ifx = -inflate_num; ifx <= inflate_num; ++ifx)
      for (int ify = -inflate_num; ify <= inflate_num; ++ify)
        for (int ifz = -0; ifz <= 0; ++ifz)
        {
          Cpt_if = Cpt + Eigen::Vector3d(ifx * _resolution, ify * _resolution,
                                         ifz * _resolution);
          Cpt_if = Cpt_if;
          racingPt.x = Cpt_if(0);
          racingPt.y = Cpt_if(1);
          racingPt.z = Cpt_if(2);
          cloudMap.push_back(racingPt);
        }
  }
}

void column_genreate_b(const Eigen::Vector3d pos, const double inflate_r)
{
  pcl::PointXYZ racingPt;
  Eigen::Vector3d Cpt;
  double inflate_num = inflate_r / _resolution;
  Cpt(0) = pos.x();
  Cpt(1) = pos.y();
  Cpt(2) = 0.0;
  for (double zm = 0; zm < pos.z(); zm += _resolution)
  {
    Cpt(2) = zm;
    Eigen::Vector3d Cpt_if;
    for (int ifx = -inflate_num; ifx <= inflate_num; ++ifx)
      for (int ify = -inflate_num; ify <= inflate_num; ++ify)
        for (int ifz = -0; ifz <= 0; ++ifz)
        {
          Cpt_if = Cpt + Eigen::Vector3d(ifx * _resolution, ify * _resolution,
                                         ifz * _resolution);
          Cpt_if = Cpt_if;
          racingPt.x = Cpt_if(0);
          racingPt.y = Cpt_if(1);
          racingPt.z = Cpt_if(2);
          cloudMap.push_back(racingPt);
        }
  }
}

void column_genreate(const Eigen::Vector3d pos1, const Eigen::Vector3d pos2, const double inflate_r, double lz, double uz)
{
  pcl::PointXYZ racingPt;
  Eigen::Vector3d Cpt;
  double inflate_num = inflate_r / _resolution;
  double l = (pos2 - pos1).norm();
  Eigen::Vector3d vec_p = (pos2 - pos1).normalized();

  Cpt(0) = pos1.x();
  Cpt(1) = pos1.y();
  Cpt(2) = 0.0;
  for (double zm = lz; zm < uz; zm += _resolution)
  {
    Cpt(0) = pos1.x();
    Cpt(1) = pos1.y();
    Cpt(2) = zm;
    double ltmp = 0;
    while (ltmp < l)
    {
      Cpt += _resolution * vec_p;
      Eigen::Vector3d Cpt_if;
      for (int ifx = -inflate_num; ifx <= inflate_num; ++ifx)
        for (int ify = -inflate_num; ify <= inflate_num; ++ify)
          for (int ifz = -0; ifz <= 0; ++ifz)
          {
            Cpt_if = Cpt + Eigen::Vector3d(ifx * _resolution, ify * _resolution,
                                           ifz * _resolution);
            racingPt.x = Cpt_if(0);
            racingPt.y = Cpt_if(1);
            racingPt.z = Cpt_if(2);
            cloudMap.push_back(racingPt);
          }
      ltmp += _resolution;
    }
  }
}

void face_genreate(const Eigen::Vector3d pos1, const Eigen::Vector3d pos2, const double inflate_r, double lz, double uz)
{
  pcl::PointXYZ racingPt;
  Eigen::Vector3d Cpt;
  double inflate_num = inflate_r / _resolution;
  double l = (pos2 - pos1).norm();
  Eigen::Vector3d vec_p = (pos2 - pos1).normalized();

  Cpt(0) = pos1.x();
  Cpt(1) = pos1.y();
  Cpt(2) = 0.0;
  for (double zm = lz; zm < uz; zm += _resolution)
  {
    Cpt(0) = pos1.x();
    Cpt(1) = pos1.y();
    Cpt(2) = zm;
    Cpt -= _resolution * vec_p;
    double ltmp = 0;
    while (ltmp < l)
    {
      Cpt += _resolution * vec_p;
      Eigen::Vector3d Cpt_if;
      for (int ifx = 10; ifx <= 10; ++ifx)
        for (int ify = -inflate_num; ify <= inflate_num; ++ify)
          for (int ifz = -0; ifz <= 0; ++ifz)
          {
            Cpt_if = Cpt + Eigen::Vector3d(ifx * _resolution, ify * _resolution,
                                           ifz * _resolution);
            racingPt.x = Cpt_if(0);
            racingPt.y = Cpt_if(1);
            racingPt.z = Cpt_if(2);
            cloudMap.push_back(racingPt);
          }
      ltmp += _resolution;
    }
  }
}

void RacingMapGenerate()
{
  pcl::PointXYZ racingPt;
  int Gatenums = Gates.cols();

  // generate circle obs
  for (int i = 0; i < Gatenums; ++i)
  {
    if (fake_gates_[i])
      continue;
    double x, y, z;
    x = Gates.col(i).x();
    y = Gates.col(i).y();
    z = Gates.col(i).z();

    double zu = z;

    x = floor(x / _resolution) * _resolution + _resolution / 2.0;
    y = floor(y / _resolution) * _resolution + _resolution / 2.0;
    z = floor(z / _resolution) * _resolution + _resolution / 2.0;
    cout << "point" << i << ": " << x << " " << y << " " << z << endl;

    double theta = gate_vectors_[i];
    Eigen::Matrix3d rotate;
    rotate << cos(theta), -sin(theta), 0.0,
        sin(theta), cos(theta), 0.0,
        0, 0, 1;

    double radius1 = inner_gate_radius_;
    double radius2 = outter_gate_radius_;

    // draw a circle centered at (x,y,z)
    int N = 128;
    double angle_step = 6.282 / double(N);

    Eigen::Vector3d cpt;
    for (double r = radius1; r < radius2; r += _resolution)
    {
      for (double angle = 0.0; angle < 6.282; angle += angle_step)
      {
        cpt(0) = 0.0;
        cpt(1) = r * cos(angle);
        cpt(2) = r * sin(angle);

        // inflate
        Eigen::Vector3d cpt_if;
        for (int ifx = -1; ifx <= 1; ++ifx)
          for (int ify = -1; ify <= 1; ++ify)
            for (int ifz = -0; ifz <= 0; ++ifz)
            {
              cpt_if = cpt + Eigen::Vector3d(ifx * _resolution, ify * _resolution,
                                             ifz * _resolution);
              cpt_if = rotate * cpt_if + Eigen::Vector3d(x, y, z);
              racingPt.x = cpt_if(0);
              racingPt.y = cpt_if(1);
              racingPt.z = cpt_if(2);
              cloudMap.push_back(racingPt);
            }
      }
    }

    double special_h = 0;
    if (i > 13)
      special_h = 4.5;
    // add 0.5m colunm
    cpt(0) = 0.0;
    cpt(1) = radius2;
    cpt(2) = 0.0;
    for (double zm = special_h; zm < zu; zm += _resolution)
    {
      cpt(2) = zm;
      Eigen::Vector3d cpt_if;
      for (int ifx = -1; ifx <= 1; ++ifx)
        for (int ify = -1; ify <= 1; ++ify)
          for (int ifz = -0; ifz <= 0; ++ifz)
          {
            cpt_if = cpt + Eigen::Vector3d(ifx * _resolution, ify * _resolution,
                                           ifz * _resolution);
            cpt_if = rotate * cpt_if + Eigen::Vector3d(x, y, 0);
            racingPt.x = cpt_if(0);
            racingPt.y = cpt_if(1);
            racingPt.z = cpt_if(2);
            cloudMap.push_back(racingPt);
          }
    }

    // add 0.5m colunm
    cpt(0) = 0.0;
    cpt(1) = -radius2;
    cpt(2) = 0.0;
    for (double z = special_h; z < zu; z += _resolution)
    {
      cpt(2) = z;
      Eigen::Vector3d cpt_if;
      for (int ifx = -1; ifx <= 1; ++ifx)
        for (int ify = -1; ify <= 1; ++ify)
          for (int ifz = -0; ifz <= 0; ++ifz)
          {
            cpt_if = cpt + Eigen::Vector3d(ifx * _resolution, ify * _resolution,
                                           ifz * _resolution);
            cpt_if = rotate * cpt_if + Eigen::Vector3d(x, y, 0);
            racingPt.x = cpt_if(0);
            racingPt.y = cpt_if(1);
            racingPt.z = cpt_if(2);
            cloudMap.push_back(racingPt);
          }
    }
  }

  Eigen::Vector3d clumn1_pos, clumn2_pos, clumn3_pos, clumn4_pos;
  clumn1_pos << 30, -13, 4;
  clumn2_pos << 23, -13, 4;
  clumn3_pos << 15, -13, 4;
  column_genreate(clumn1_pos, 0.75);
  column_genreate(clumn2_pos, 0.75);
  column_genreate(clumn3_pos, 0.75);

  double cps[] = {3, -5, -20, -28};

  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j <= 8; j++)
    {
      Eigen::Vector3d tmp_p(j * 5, cps[i], 4.5);
      column_genreate_b(tmp_p, 0.5);
    }
  }

  for (int i = 0; i < 2; i++)
  {
    for (int j = 1; j <= 8; j++)
    {
      Eigen::Vector3d tmp_p1((j - 1) * 5, cps[i], 4.5);
      Eigen::Vector3d tmp_p2(j * 5, cps[i], 4.5);
      column_genreate(tmp_p1, tmp_p2, 0.5, 4.3, 4.5);
    }
  }

  Eigen::Vector3d tmp_p1(-1, -24, 4.5);
  Eigen::Vector3d tmp_p2(39.5, -24, 4.5);
  face_genreate(tmp_p1, tmp_p2, 4.5, 4.3, 4.5);

  cloudMap.width = cloudMap.points.size();
  cloudMap.height = 1;
  cloudMap.is_dense = true;

  _map_ok = true;

  ROS_WARN("Finished generate random map ");
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "random_map_sensing");
  ros::NodeHandle n("~");

  _all_map_pub = n.advertise<sensor_msgs::PointCloud2>("global_cloud", 1);
  // map_inflate_pub_ = n.advertise<sensor_msgs::PointCloud2>("occupancy_inflate", 10);
  circle_pub_ = n.advertise<geometry_msgs::PoseStamped>("/circle_odom", 1);
  hPolyPub_ = n.advertise<decomp_ros_msgs::PolyhedronArray>("polyhedra", 1);
  global_traj_pub_ = n.advertise<quadrotor_msgs::PolynomialTrajectory>("global_trajectory", 1);
  gridmap_vs_pub = n.advertise<sensor_msgs::PointCloud2>("vs_gridmap", 1);
  gridmap_inflate_vs_pub = n.advertise<sensor_msgs::PointCloud2>("vs_gridmap_inflate", 1);

  click_map_pub_ = n.advertise<sensor_msgs::PointCloud2>("/pcl_render_node/local_map", 1);
  // ros::Subscriber click_sub = n.subscribe("/goal", 10, clickCallback);
  n.param("inner_gate_radius", inner_gate_radius_, 1.0);
  n.param("outter_gate_radius", outter_gate_radius_, 1.0);

  n.param("init_state_x", _init_x, 0.0);
  n.param("init_state_y", _init_y, 0.0);
  n.param("init_state_z", _init_z, 1.0);

  n.param("map/x_size", _x_size, 50.0);
  n.param("map/y_size", _y_size, 50.0);
  n.param("map/z_size", _z_size, 5.0);
  n.param("map/resolution", _resolution, 0.1);
  n.param("map/inflate_size", inflate_size_, 0.3);

  n.param("ObstacleShape/theta", theta_, 7.0);

  n.param("sensing/radius", _sensing_range, 10.0);
  n.param("sensing/rate", _sense_rate, 10.0);

  n.param("corridor/l_corridor", l_corridor_, 2.0);
  n.param("corridor/bbox_width", bbox_width_, 4.0);
  n.param("corridor/desired_len", desired_len_, 2.0);

  n.getParam("gate_oritation", gate_vectors_);

  std::vector<double> tmp_gate_pos_vec;

  n.getParam("gate_pos_vec", tmp_gate_pos_vec);

  std::cout << "tmp_gate_pos_vec.size(): " << tmp_gate_pos_vec.size() << std::endl;
  int num = tmp_gate_pos_vec.size() / 3;
  std::cout << "num: " << num << std::endl;

  n.getParam("target_pos", target_p_vec_);
  target_p_ << target_p_vec_[0], target_p_vec_[1], target_p_vec_[2];

  Gates.resize(3, num);
  for (int i = 0; i < num; i++)
  {
    Gates.col(i) << tmp_gate_pos_vec[3 * i], tmp_gate_pos_vec[3 * i + 1], tmp_gate_pos_vec[3 * i + 2];
  }

  std::vector<double> special_gate_index;
  n.getParam("special_gate_index", special_gate_index);

  for (int i = 0; i < num; i++)
  {
    bool is_fake_gate = false;
    for (int j = 0; j < int(special_gate_index.size()); j++)
    {
      if (special_gate_index[j] == i)
      {
        is_fake_gate = true;
        break;
      }
    }
    fake_gates_.push_back(is_fake_gate);
  }

  ROS_WARN("Initializae map");

  visPtr_ = std::make_shared<visualization::Visualization>(n);

  ros::Duration(0.5).sleep();
  ROS_WARN("Start generate map");

  RacingMapGenerate();

  ROS_WARN("finish generate map");

  ros::Rate loop_rate(_sense_rate);

  ROS_WARN("start publish map");

  pcl::toROSMsg(cloudMap, globalMap_pcd);
  globalMap_pcd.header.frame_id = "world";
  wr_msg::writeMsg(globalMap_pcd, ros::package::getPath("plan_manage") + "/offlines/globalmap.bin");

  pcl::io::savePCDFileBinary(ros::package::getPath("plan_manage") + "/pcd/globalmap.pcd", cloudMap);

  ROS_WARN("finish publish map");

  while (ros::ok())
  {
    _all_map_pub.publish(globalMap_pcd);
    ros::spinOnce();
    loop_rate.sleep();
  }
}