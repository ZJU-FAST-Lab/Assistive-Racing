#pragma once

#include<Eigen/Dense>
#include<Eigen/Geometry>
#include<iostream>
#include<ros/ros.h>
#include "uav_utils/geometry_utils.h"
// #include "trajectory_utils/poly_traj_utils.hpp"

using namespace std;

struct Action{
  double vx;
  double vz;
  double omega;
};

namespace RC{
class RC_Op {
  ros::NodeHandle nh_;
  Eigen::VectorXi channels;
  Action input_act;
  Eigen::Vector3d joy_acc_cmd;
  Eigen::Vector3d joy_vel_cmd;  // [vel_x, vel_y, vel_z, yaw_dot]
  Eigen::Quaterniond q_odom;
  // Eigen::MatrixXd q_state;
  // Trajectory primitive;

  double PI = 3.1415926;
  double R,P,Y,thrust;
  double Y_dot_max_, R_max_, P_max_, thr_max_;
  double Y_cmd_dot_max_;
  double rc_velx_max_, rc_vely_max_, rc_velz_max_;

  double Cmd_Ydot;

  bool if_Auto_Yaw_;

  
public:
  RC_Op(ros::NodeHandle& nh);
  ~RC_Op(){};

  void CalculateAccCmd();

  void CalculateVelCmd();

  void Inputs(Eigen::VectorXi chan);

  Eigen::Vector3d GetJoyAccCmd();

  Eigen::Vector3d GetJoyVelCmd();

  double GetYaw();

  double GetYawCmd();
  
  void SetYaw(const double &Yaw); 
};

};


