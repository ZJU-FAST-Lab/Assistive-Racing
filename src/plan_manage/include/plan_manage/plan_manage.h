#ifndef _PLAN_MANAGE_H_
#define _PLAN_MANAGE_H_

#include "ros/ros.h"
#include <ros/package.h>
#include <math.h>
#include <iostream>
#include <Eigen/Eigen>

#include <nav_msgs/Odometry.h>
#include <std_msgs/Float32.h>
#include <geometry_msgs/PoseStamped.h>
#include <quadrotor_msgs/DebugData.h>
#include <quadrotor_msgs/RiskDebug.h>
#include <quadrotor_msgs/PolynomialTrajectory.h>
#include <quadrotor_msgs/HpolyTimeSerial.h>
#include <geometry_msgs/Point.h>
#include <sensor_msgs/Joy.h>
#include <mavros_msgs/RCIn.h>
#include <sensor_msgs/Range.h>
#include <tf2_ros/transform_broadcaster.h>
#include <tf2_geometry_msgs/tf2_geometry_msgs.h>

#include "traj_opt/lc_yaw_traj_opt_varied.h"
#include "traj_opt/uniform_minco_progress.h"
#include "traj_opt/uniform_minco_emergency.h"

#include "RC/RC.h"
#include "risk_opt/risk_monitor.h"
#include "corridor_helper/corridor_helper.h"

#include <decomp_ros_utils/data_ros_utils.h>
#include <decomp_util/ellipsoid_decomp.h>

#include "wr_msg/wr_msg.hpp"

namespace planner
{

  enum Traj_STATE
  {
    PRIM,
    ASSIST
  };

  class PlanManager
  {
  private:
    ros::NodeHandle nh_;
    //------------------------------------------------------------------- ros related -------------------------------------------------------------------
    ros::Subscriber odom_sub_, circle_sub_, depth_img_sub_, RC_sub_, polyhedra_sub_, global_traj_sub_;
    ros::Publisher traj_pub_, desired_yaw_pub_, hPolyPub_, OptTime_pub_, Debug_pub, RiskDebug_pub;
    ros::Timer process_timer_, monitor_timer_, checker_timer_;
    ros::Time start_timer;

    // -------------------------------------------------------------------state related -------------------------------------------------------------------
    Eigen::Vector3d odom_p_, odom_v_;
    Eigen::Vector3d start_p_, start_v_, start_a_, start_j_;
    Eigen::Vector3d target_p_, target_v_, target_a_, target_j_;
    Eigen::Vector3d global_target_p_, global_yaw_target_p_;

    //------------------------------------------------------------------- replan related -------------------------------------------------------------------
    double min_replan_dist_, min_replan_duration_, min_circle_update_dist_;
    double replan_duration_, alpha_gate_update_;
    double replan_horizon_;

    bool get_ref_traj_ = false;

    int ref_replan_index_;

    double init_x_, init_y_, init_z_;
    double init_yaw_;

    int current_id_;

    bool fly_start_trgger_ = false;

    std::vector<std::pair<double, double>> time_serial_;

    int N_1_;

    //-------------------------------------------------------------------trajectory related -------------------------------------------------------------------
    int gate_index_ = -1;
    bool get_glb_traj_ = false, get_lc_traj_ = false, get_lc_prim_traj_ = false;
    bool get_lc_yaw_traj_ = false, get_glb_yaw_traj_ = false;
    Trajectory glb_traj_, traj_emergency_;
    Trajectory7 lc_traj_, lc_prim_traj_;
    Trajectory1D glb_yaw_traj_, lc_yaw_traj_;
    Trajectory1D lc_time_traj_;
    ros::Time glb_start_timestamp_, lc_start_timestamp_, lc_prim_start_timestamp_, lc_ref_start_timestamp_, replan_timestamp_emergency_;
    ros::Time lc_yaw_start_timestamp_, glb_yaw_start_timestamp_;

    // -------------------------------------------------------------------gete related -------------------------------------------------------------------
    std::vector<Eigen::Vector3d> gate_pos_list_;
    std::vector<Eigen::Vector3d> gate_oritation_vector;
    double gate_radius_;

    // gate_pos
    Eigen::Vector3d circle_pos_;
    bool circle_pos_receive_ = false;

    // ------------------------------------------------------------------- Ptrs -------------------------------------------------------------------
    traj_opt::Config config_;
    std::shared_ptr<visualization::Visualization> visPtr_;

    std::shared_ptr<traj_opt::LocalYawOpt> lcYawTrajOptPtr_;
    std::shared_ptr<traj_opt::UniformMincoProgress> assistProgressOptPtr_;
    std::shared_ptr<traj_opt::UniformMincoStop> EmergencyStopOptPtr_;
    std::shared_ptr<corridor_opt::CorridorHelper> CorridorPtr_;
    std::shared_ptr<RC::RC_Op> RcPtr_;
    std::shared_ptr<Risk_opt::risk_monitor> RiskPtr_;

    //-------------------------------------------------------------------RC_related -------------------------------------------------------------------
    Trajectory7 RC_prim;
    Trajectory1D Yaw_prim;
    Traj_STATE traj_state_, last_traj_state_;
    bool ifAutochange_ = false;

    double rc_expect_t_;
    bool if_Auto_Yaw_ = false;

    //-------------------------------------------------------------------risk related -------------------------------------------------------------------
    ros::Time t_risk_trigger_, t_prim_trigger_;
    double progress_deviative_up_threshold_, progress_deviative_low_threshold_;
    double rc_deviative_up_threshold_, rc_deviative_low_threshold_;

    // -------------------------------------------------------------------aux realted  -------------------------------------------------------------------
    double tmax;

    // new hpolys recieved
    bool receive_hpolys_ = false;
    std::vector<Eigen::MatrixXd> TunnelLikeHpolys_;

    // Debug data
    quadrotor_msgs::DebugData Debug_msg_;
    quadrotor_msgs::RiskDebug Risk_msg_;

    // to determine the finish
    Eigen::Vector3d terminal_pos_, terminal_vec_;
    bool reach_flag_ = false;
    bool if_emergency_stop_ = false;

    // msg reader
    std::string glb_traj_path_, corridor_path_, hpoly_timeserial_path_;

  public:
    PlanManager(ros::NodeHandle &nh);
    ~PlanManager() {}

    double getYawVecToGate(const Eigen::Vector3d pos, const int current_ids);

    void ReadGlobalTraj();

    void ReadCorridors();

    void SetEmergencyStopBox();

    void StateChange(Traj_STATE state);

    void PrintState(Traj_STATE state);

    void SafeTrajCallback(const ros::TimerEvent &event);

    void TicCallback(const ros::TimerEvent &event);

    void polyhedra_callback(const decomp_ros_msgs::PolyhedronArrayPtr &msgPtr);

    void GlobalTraj_callback(const quadrotor_msgs::PolynomialTrajectoryPtr &msgPtr);

    void rcvOdometryCallback(const nav_msgs::Odometry &odom);

    void rcvCircleCallback(const geometry_msgs::PoseStamped &msg);

    void rcvJoyCmdCallback(const mavros_msgs::RCIn::ConstPtr &msg);

    void local_replan(const int &gate_index);

    void get_replan_state_prim(ros::Time &replan_timestamp, Eigen::MatrixXd &initState, Trajectory7 &tmp_traj, double &replan_t);

    void get_yaw_replan_state(ros::Time &replan_timestamp, Eigen::VectorXd &initState, Eigen::VectorXd &finState, Trajectory1D &tmp_traj, const int &gate_index, double &replan_t);

    Eigen::MatrixXd StateIntegral(const Eigen::MatrixXd State, const double integralT);

    void get_replan_state(ros::Time &replan_timestamp, Eigen::MatrixXd &initState, Eigen::MatrixXd &finState, Trajectory7 &tmp_traj, const int &gate_index, double &replan_t);

    void get_replan_state(ros::Time &replan_timestamp, Eigen::MatrixXd &initState, Eigen::MatrixXd &finState, Trajectory &tmp_traj, const int &gate_index, double &replan_t);

    void visualize_traj();

    void publish_traj(const Trajectory &traj, const ros::Time &traj_timestamp);

    void publish_traj(const Trajectory7 &traj, const ros::Time &traj_timestamp);

    quadrotor_msgs::PolynomialTrajectory traj2msg(const Trajectory &traj, const ros::Time &traj_timestamp);

    quadrotor_msgs::PolynomialTrajectory traj2msg(const Trajectory7 &traj, const ros::Time &traj_timestamp);

    void publish_desired_yaw();

    Trajectory calculatePrimitive(Eigen::Vector3d acc_cmd,
                                  Eigen::Vector3d q_pos,
                                  Eigen::Vector3d q_vel,
                                  double dur,
                                  bool isAcc = true);

    Trajectory7 calculateVelBasedPrimitive(Eigen::Vector3d v_cmd,
                                           Eigen::Vector3d q_pos,
                                           Eigen::Vector3d q_vel,
                                           Eigen::Vector3d q_acc,
                                           Eigen::Vector3d q_jrk,
                                           double dur);

    Trajectory calculateStopPrimitive(Eigen::Vector3d end_pos,
                                      Eigen::Vector3d q_pos,
                                      Eigen::Vector3d q_vel,
                                      Eigen::Vector3d q_acc,
                                      double dur);

    Trajectory1D calculateYawPrimitive(const Eigen::Vector3d yaw_init,
                                       const double yaw_dot_cmd,
                                       const double dur);

    inline void visCorridor(const vec_E<Polyhedron3D> &polyhedra);
    inline void visCorridor(const std::vector<Eigen::MatrixXd> &hPolys);

  }; // class plan_manager

} // namespace planner

#endif
