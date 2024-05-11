#pragma once

#include<Eigen/Dense>
#include<Eigen/Geometry>
#include<ros/ros.h>
#include "uav_utils/geometry_utils.h"
#include "trajectory_utils/poly_traj_utils.hpp"
#include <optimization_opt/lbfgs_raw.hpp>
#include "corridor_helper/corridor_helper.h"

using namespace std;

namespace Risk_opt{
class risk_monitor{
  ros::NodeHandle nh_;
  int N_, K_;
  bool lock_;
  double beta;
  Trajectory ref_traj, evaluate_traj;
  ros::Time ref_start_time_, eval_start_time_;
  bool if_ref_receive_, if_eval_receive_;
  bool if_state_receive_;

  double deviate_cost_;
  double rc_deviate_cost_;
  double time_cost_;
  double acc_cost_;

  double time_threshold_;
  double deviate_thrushold_;
  double acc_threshold_;

  double threshold_for_assist_, threshold_for_prim_;

  Eigen::Vector3d next_gate_pos_;

  Eigen::VectorXd x_;
  double risk_vel_max_, risk_acc_max_;   
  double BigPen_; 
  Eigen::MatrixXd current_state_;

  Eigen::MatrixXd hpoly1_, hpoly2_;

  lbfgs::lbfgs_parameter_t lbfgs_params_;

  Eigen::MatrixXd b_;

  double minSafeT_;
  double tol_;

  double rhoT_;

  double current_id_;

  std::shared_ptr<corridor_opt::CorridorHelper> CorridorPtr_;

public:
  risk_monitor(ros::NodeHandle nh){
      nh_ = nh;
      lock_ = false;
      nh_.getParam("N_1", N_);
      nh_.getParam("time_threshold", time_threshold_);
      nh_.getParam("deviate_thrushold", deviate_thrushold_);
      nh_.getParam("acc_threshold", acc_threshold_);
      nh_.getParam("threshold_for_assist", threshold_for_assist_);
      nh_.getParam("threshold_for_prim", threshold_for_prim_);

      nh_.param<double>("risk_vel_max", risk_vel_max_, 7.0);
      nh_.param<double>("risk_acc_max", risk_acc_max_, 7.0);

      nh_.param<double>("minSafeT", minSafeT_, 0.3);

      if_ref_receive_ = false;
      if_eval_receive_ = false;
      if_state_receive_ = false;
      current_state_.resize(3, 4);

      x_.resize(1);
      lbfgs_params_.mem_size = 16;
      lbfgs_params_.past = 3;
      lbfgs_params_.g_epsilon = 0.0;
      lbfgs_params_.min_step = 1e-32;
      lbfgs_params_.delta = 1.0e-4;

      b_.resize(6,3);
      b_.setZero();

      BigPen_ = 100000.0;
      tol_ = 1e-8;

      K_ = 5;
      rhoT_ = 10000.0;
      deviate_cost_ = 0;

      cout<<"--------------------Risk Monitor Params-------------------------"<<endl;
      cout<<" N_1: "<<N_<<" time_threshold_: "<<time_threshold_<<" deviate_thrushold_: "<<deviate_thrushold_<<" acc_thrushold_: "<<acc_threshold_<<endl;
      cout<<"threshold_for_assist_: "<<threshold_for_assist_<<" threshold_for_prim_: "<<threshold_for_prim_<<endl;
      cout<<"risk_vel_max_: "<<risk_vel_max_<<" risk_acc_max_: "<<risk_acc_max_<<" minSafeT: "<< minSafeT_ <<endl;
  };
  ~risk_monitor(){
      // delete x_;
  };

  void setRefTraj(Trajectory& traj, ros::Time t){
      ref_traj = traj;
      ref_start_time_ = t;
      if_ref_receive_ = true;
  }

  void setEvaluateTraj(Trajectory& traj, ros::Time t){
      evaluate_traj = traj;
      eval_start_time_ = t;
      if_eval_receive_ = true;
  }

  void setDeviateCost(double cost){
      deviate_cost_ = cost;
  }

  void setRcDeviateCost(double cost){
      rc_deviate_cost_ = cost;
  }

  double getDeviateCost(){
      return deviate_cost_;
  }

  double getRcDeviateCost(){
      return rc_deviate_cost_;
  }   

  void setCurrentState(const Eigen::MatrixXd state){
      current_state_ = state; 
      b_.block(0,0,3,3) = state.block(0, 0, 3, 3).transpose();
      if_state_receive_ = true;
  }

  void setGatePose(Eigen::Vector3d next_gate_pos){
      next_gate_pos_ = next_gate_pos;
  }

  void setCorridorHelper(const std::shared_ptr<corridor_opt::CorridorHelper>& CorridorPtr) {
      CorridorPtr_ = CorridorPtr;
  }

  void setCurrentCorridorId(const int idx){
      current_id_ = idx;
  }

  double getDeviateCost(ros::Time eval_t);

  double time_evaluate(ros::Time eval_t);    

  double getProgress(Eigen::Vector3d pos);

  bool VelvecCheck(double cos_threshold_);

  double SeekMinVelvec();

  // bool ifAssist(){
  //     if(time_cost_ > time_threshold_ || acc_cost_ > acc_threshold_) return true;
  //     return false;
  // }

  bool ifAssist();
  bool ifUserControl();

  bool SafetyCheck(Trajectory &emergency_traj);

  static double objectiveFunc(void *ptrObj, const Eigen::VectorXd &x, Eigen::VectorXd& grad);

  bool MaxConstraint(const Eigen::Vector3d &derivative,
                      Eigen::Vector3d &grad,
                      double &cost,
                      const double &max_value);

  bool ifCollision(const Trajectory &traj, const CoefficientMat &coeff, const double &duration);

};
}   // namespace Risk_opt