#include "RC/RC.h"
using namespace std;

namespace RC
{

  RC_Op::RC_Op(ros::NodeHandle &nh)
  {
    int channel_size_;
    nh_ = nh;
    nh_.param<int>("channel_size", channel_size_, 8.0);
    channels.resize(channel_size_);

    nh_.param<double>("R_max", R_max_, 0.5);
    nh_.param<double>("P_max", P_max_, 0.5);
    nh_.param<double>("Y_dot_max", Y_dot_max_, 0.005);
    nh_.param<double>("Y_cmd_dot_max", Y_cmd_dot_max_, 0.5);
    nh_.param<double>("RCthracc_max", thr_max_, 32.0);
    nh_.param<double>("rc_velx_max", rc_velx_max_, 7.0);
    nh_.param<double>("rc_vely_max", rc_vely_max_, 5.0);
    nh_.param<double>("rc_velz_max", rc_velz_max_, 3.0);
    nh_.param<bool>("if_Auto_Yaw", if_Auto_Yaw_, true);

    R_max_ = PI * R_max_;
    P_max_ = PI * P_max_;
    Y_dot_max_ = PI * Y_dot_max_;
    Y_cmd_dot_max_ = PI * Y_cmd_dot_max_;

    Y = 0;
    R = 0;
    P = 0;
    thrust = 0;

    joy_acc_cmd.setZero();
    joy_vel_cmd.setZero();

    cout << "--------------------RC Params-------------------------" << endl;
    cout << " R_max: " << R_max_ << " Y_dot_max: " << Y_dot_max_ << " P_max: " << P_max_ << " thr_max: " << thr_max_ << endl
         << "Y_cmd_dot_max: " << Y_cmd_dot_max_ << std::endl;
    cout << "rc_velx_max: " << rc_velx_max_ << " rc_vely_max: " << rc_vely_max_ << " rc_velz_max: " << rc_velz_max_ << endl;
    cout << "if_Auto_Yaw" << if_Auto_Yaw_ <<std::endl;
  }

  void RC_Op::Inputs(Eigen::VectorXi chan)
  {
    assert(channels.size() == chan.size());
    channels = chan;
    CalculateVelCmd();
  }

  Eigen::Vector3d RC_Op::GetJoyAccCmd()
  {
    return joy_acc_cmd;
  }

  Eigen::Vector3d RC_Op::GetJoyVelCmd()
  {
    return joy_vel_cmd;
  }

  double RC_Op::GetYaw()
  {
    return Y;
  }

  void RC_Op::SetYaw(const double &Yaw)
  {
    Y = Yaw;
  }

  double RC_Op::GetYawCmd()
  {
    return Cmd_Ydot;
  }

  void RC_Op::CalculateAccCmd()
  {
    // work well
    // Eigen::Vector3d odom_euler = uav_utils::R_to_ypr(q_odom.toRotationMatrix());
    // double odom_yaw = odom_euler(0);
    R = channels(0) - 1500;
    P = channels(1) - 1500;
    thrust = channels(2) - 1000;

    // caculate the R from B to W
    Eigen::Vector3d eulerAngle;
    eulerAngle(0) = R_max_ * R / 500;
    eulerAngle(1) = P_max_ * P / 500;
    eulerAngle(2) = Y;

    Eigen::AngleAxisd rollAngle(Eigen::AngleAxisd(eulerAngle(0), Eigen::Vector3d::UnitX()));
    Eigen::AngleAxisd pitchAngle(Eigen::AngleAxisd(eulerAngle(1), Eigen::Vector3d::UnitY()));
    Eigen::AngleAxisd yawAngle(Eigen::AngleAxisd(eulerAngle(2), Eigen::Vector3d::UnitZ()));

    Eigen::Quaterniond q;
    q = yawAngle * pitchAngle * rollAngle;

    joy_acc_cmd = q.toRotationMatrix() * Eigen::Vector3d(0, 0, thrust / 1000 * thr_max_);
  }

  void RC_Op::CalculateVelCmd()
  {
    double Y_dot = double(channels(3)) - 1500;

    double forward_x = (double(channels(1)) - 1500) / 500 * rc_velx_max_;
    double forward_y = -(double(channels(0) - 1500)) / 500 * rc_vely_max_;
    double upward = (double(channels(2) - 1500)) / 500 * rc_velz_max_;

    Eigen::Vector3d v_cmd;
    v_cmd << forward_x, forward_y, upward;

    if (!if_Auto_Yaw_)
      Y = Y + Y_cmd_dot_max_ * Y_dot / 1000;
    Cmd_Ydot = Y_dot_max_ * Y_dot / 1000;

    Eigen::AngleAxisd yawAngle(Eigen::AngleAxisd(Y, Eigen::Vector3d::UnitZ()));
    Eigen::Quaterniond q(yawAngle);

    joy_vel_cmd = q.toRotationMatrix() * v_cmd;
  }
};
