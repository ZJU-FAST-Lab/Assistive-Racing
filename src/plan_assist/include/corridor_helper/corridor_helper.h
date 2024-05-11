#pragma once

#include<Eigen/Dense>
#include<Eigen/Geometry>
#include<iostream>
#include<ros/ros.h>
// #include "uav_utils/geometry_utils.h"
#include "trajectory_utils/poly_traj_utils.hpp"

using namespace std;

namespace corridor_opt{
class CorridorHelper {
  ros::NodeHandle nh_;
  std::vector<Eigen::MatrixXd> corridors_;
  int num_corridor;
  bool if_hpolys_receive_;
  
public:
  CorridorHelper(ros::NodeHandle& nh){
    nh_ = nh;
    if_hpolys_receive_ = false;
  };
  ~CorridorHelper(){};

  void setCorridor(const std::vector<Eigen::MatrixXd> corridors){
    corridors_ = corridors;
    num_corridor = corridors_.size();
    if_hpolys_receive_ = true;
  }

  int getCorridorNum(){
    return num_corridor;
  }

  Eigen::MatrixXd getHpoly(const int idx){
    return corridors_[idx];
  }

  bool ifCorridorReceived(){
    return if_hpolys_receive_;
  }

  bool isInhpolys(const Eigen::Vector3d pos, const Eigen::MatrixXd hpoly, double &min_d);

  bool isInhpolysConversative(const Eigen::Vector3d pos, const Eigen::MatrixXd hpoly, double &min_d);

  int getBestHpoly(const Eigen::Vector3d pos, const int huerist_idx);

  int getBestHpoly(const Eigen::Vector3d pos, const int huerist_idx, const int expected_forward_num);

  double getBestDistInHpoly(const Eigen::Vector3d pos, const int huerist_idx);

  double getBestDistInHpoly(const Eigen::Vector3d pos, const int huerist_idx, const int expected_forward_num);

  int getCurrentHpolyForward(const Eigen::Vector3d pos, const int huerist_idx);

  int getCurrentHpolyBackward(const Eigen::Vector3d pos, const int huerist_idx);

  int SafeCheckHpolyForward(const Eigen::Vector3d pos, const int huerist_idx);

  int SafeCheckHpolyBackward(const Eigen::Vector3d pos, const int huerist_idx);

  int getfarthestHpolyForward(const Eigen::Vector3d pos, const int huerist_idx, const int expected_forward_num);

  int SafeCheckHpolyRange(const Eigen::Vector3d pos, const int huerist_idx, int lower, int upper);

};

};


