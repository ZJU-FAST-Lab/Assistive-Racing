#include <iostream>
#include <ros/ros.h>
#include <Eigen/Eigen>

#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/Vector3.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/PointCloud2.h>

#include <ros/package.h>
#include "wr_msg/wr_msg.hpp"

ros::Publisher circle_pub_, glbpcd_pub_;
ros::Subscriber odom_sub_;

bool has_odom_= false;
int idx_now_ = 0;
double gate_radius_;
Eigen::Vector3d state_pos_;
Eigen::Vector3d target_p_;

std::vector<Eigen::Vector3d> gate_pos_list_;
std::vector<Eigen::Vector3d> gate_vec_list_;

std::vector<bool> fake_gates_;
std::string map_path_;

void pubCurrentGate(){
  if(!has_odom_) return;
  geometry_msgs::PoseStamped gate_pos_msg;
  if(idx_now_ >= int(gate_pos_list_.size())){
    // ROS_WARN("Pass through all gates");
    return;
  } 
  Eigen::Vector3d tmp = gate_pos_list_[idx_now_];
  Eigen::Vector3d theta_tmp = gate_vec_list_[idx_now_];
  double dist1 = (state_pos_-tmp).norm();
  double dist = (tmp - state_pos_).dot(theta_tmp);

  if(fake_gates_[idx_now_]){
    if(dist < -1e-6){
      idx_now_++;
    } 
  }
  else if(dist < -1e-6 && dist1 < gate_radius_){
     idx_now_++;
  }
  gate_pos_msg.pose.position.x = tmp.x();
  gate_pos_msg.pose.position.y = tmp.y();
  gate_pos_msg.pose.position.z = tmp.z();
  circle_pub_.publish(gate_pos_msg);
}

void rcvOdometryCallbck(const nav_msgs::Odometry odom) {
  if (odom.child_frame_id == "X" || odom.child_frame_id == "O") return;
  has_odom_ = true;

  state_pos_ << odom.pose.pose.position.x, odom.pose.pose.position.y, odom.pose.pose.position.z;
}

int main(int argc, char** argv) {
  ros::init(argc, argv, "fake_circle_node");
  ros::NodeHandle nh("~");
  
  circle_pub_ = nh.advertise<geometry_msgs::PoseStamped>("circle_odom",1 );
  glbpcd_pub_ = nh.advertise<sensor_msgs::PointCloud2>("vis_gridmap", 1);

  odom_sub_ = nh.subscribe("odometry", 50, rcvOdometryCallbck);

  std::vector<double> tmp_gate_pos_vec;
  std::vector<double> tmp_gate_oritation_vec;
  std::vector<double> tmp_target_p_vec;
  std::vector<double> special_gate_index;

  nh.getParam("gate_oritation", tmp_gate_oritation_vec);
  nh.getParam("gate_pos_vec", tmp_gate_pos_vec);
  nh.getParam("target_pos", tmp_target_p_vec);
  nh.getParam("special_gate_index", special_gate_index);
  nh.param("gate_radius", gate_radius_, 1.0);
  nh.param<std::string>("map_path", map_path_, std::string(ros::package::getPath("plan_manage") + "/offlines/globalmap.bin"));

  
  int num = tmp_gate_pos_vec.size() / 3;
  std::cout<<"[fake_circle_node] num: "<<num<<std::endl;
  
  target_p_ << tmp_target_p_vec[0], tmp_target_p_vec[1], tmp_target_p_vec[2];

  gate_pos_list_.resize(num);
  gate_vec_list_.resize(num);

  for(int i = 0; i < num; i++) {
    gate_pos_list_[i] << tmp_gate_pos_vec[3 * i], tmp_gate_pos_vec[3 * i + 1], tmp_gate_pos_vec[3 * i + 2];
    gate_vec_list_[i] << cos(tmp_gate_oritation_vec[i]), sin(tmp_gate_oritation_vec[i]), 0;
  }

  for(int i = 0; i< num; i++){
    bool is_fake_gate = false;
    for(int j = 0; j < int(special_gate_index.size()); j++){
      if(special_gate_index[j] == i){
        is_fake_gate = true;
        break;
      }
    }
    fake_gates_.push_back(is_fake_gate);
  }

  sensor_msgs::PointCloud2 globalMap_pcd;
  ROS_WARN("[fake_circle_node] Reading the map path from %s", map_path_.c_str());
  wr_msg::readMsg(globalMap_pcd, map_path_);

  
  ROS_WARN("[fake_circle_node] start pub current gate!");
  ros::Rate loop_rate(100);
  static ros::Time t1 = ros::Time::now();
  int cnt = 0;
  while (ros::ok()) {

    pubCurrentGate();

    // only pub the global map pcd for 10 times
    if((ros::Time::now() - t1).toSec() > 1.0 && cnt < 10){
      glbpcd_pub_.publish(globalMap_pcd);
      t1 = ros::Time::now();
      cnt ++;
    }

    ros::spinOnce();
    loop_rate.sleep();
  }

  return 0;
}