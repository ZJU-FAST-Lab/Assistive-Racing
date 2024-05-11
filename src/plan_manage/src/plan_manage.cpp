#include <plan_manage/plan_manage.h>

namespace planner
{
  // auxilury functions
  static double sign(double a)
  {
    if (a > 0)
      return 1;
    else
      return -1;
  }

  static double find_angle(double a, double b, double &sign1)
  {
    if (a > b)
    {
      while ((a - b) > 3.1415926)
      {
        b += 2 * 3.1415926;
      }
      sign1 = sign(a - b);
      return abs(a - b);
    }
    else
    {
      while ((b - a) > 3.1415926)
      {
        a += 2 * 3.1415926;
      }
      sign1 = sign(a - b);
      return abs(a - b);
    }
  }

  static void roundAngle(double &a, double &b)
  {
    if (abs(a - b) > 3.14)
    {
      double signal;
      double a1 = find_angle(a, b, signal);
      a = b + signal * a1;
    }
    return;
  }

  // visualization realted functions
  void PlanManager::visCorridor(const vec_E<Polyhedron3D> &polyhedra)
  {
    decomp_ros_msgs::PolyhedronArray poly_msg = DecompROS::polyhedron_array_to_ros(polyhedra);
    poly_msg.header.frame_id = "world";
    poly_msg.header.stamp = ros::Time::now();
    hPolyPub_.publish(poly_msg);
  }

  void PlanManager::visCorridor(const std::vector<Eigen::MatrixXd> &hPolys)
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

  PlanManager::PlanManager(ros::NodeHandle &nh) : nh_(nh)
  {
    circle_sub_ = nh_.subscribe("circle_odom", 1, &PlanManager::rcvCircleCallback, this);
    odom_sub_ = nh_.subscribe("odom", 1, &PlanManager::rcvOdometryCallback, this);
    RC_sub_ = nh_.subscribe("joycmd", 1, &PlanManager::rcvJoyCmdCallback, this);

    traj_pub_ = nh_.advertise<quadrotor_msgs::PolynomialTrajectory>("trajectory", 50);
    desired_yaw_pub_ = nh_.advertise<std_msgs::Float32>("desired_yaw", 1);
    OptTime_pub_ = nh_.advertise<std_msgs::Float32>("OptTime", 1);
    Debug_pub = nh_.advertise<quadrotor_msgs::DebugData>("DebugData", 1);
    RiskDebug_pub = nh_.advertise<quadrotor_msgs::RiskDebug>("RiskDebugData", 1);
    hPolyPub_ = nh_.advertise<decomp_ros_msgs::PolyhedronArray>("polyhedra", 100);

    config_.load(nh_);
    process_timer_ = nh_.createTimer(ros::Duration(config_.safe_traj_gen_time_), &PlanManager::SafeTrajCallback, this);
    monitor_timer_ = nh_.createTimer(ros::Duration(0.01), &PlanManager::TicCallback, this);

    // shared ptr
    CorridorPtr_ = std::make_shared<corridor_opt::CorridorHelper>(nh_);
    visPtr_ = std::make_shared<visualization::Visualization>(nh_);

    lcYawTrajOptPtr_ = std::make_shared<traj_opt::LocalYawOpt>(nh_, config_);
    lcYawTrajOptPtr_->setVisualizer(visPtr_);

    assistProgressOptPtr_ = std::make_shared<traj_opt::UniformMincoProgress>(nh_, config_);
    assistProgressOptPtr_->setVisualizer(visPtr_);
    assistProgressOptPtr_->setCorridorHelper(CorridorPtr_);

    EmergencyStopOptPtr_ = std::make_shared<traj_opt::UniformMincoStop>(nh_, config_);
    EmergencyStopOptPtr_->setVisualizer(visPtr_);

    RcPtr_ = std::make_shared<RC::RC_Op>(nh_);
    RiskPtr_ = std::make_shared<Risk_opt::risk_monitor>(nh_);
    RiskPtr_->setCorridorHelper(CorridorPtr_);

    visPtr_->registe<nav_msgs::Path>("global_traj");
    visPtr_->registe<nav_msgs::Path>("local_traj");
    visPtr_->registe<sensor_msgs::PointCloud2>("global_traj_wayPts");
    visPtr_->registe<sensor_msgs::PointCloud2>("local_traj_wayPts");

    nh_.param<double>("replan_duration", replan_duration_, 0.1);
    nh_.param<double>("replan_horizon", replan_horizon_, 5.0);
    nh_.param<double>("alpha_gate_update", alpha_gate_update_, 1.0);
    nh_.param<double>("min_replan_dist", min_replan_dist_, 1.0);
    nh_.param<double>("min_replan_duration", min_replan_duration_, 0.15);
    nh_.param<double>("min_circle_update_dist", min_circle_update_dist_, 2.0);
    nh_.param<double>("progress_deviative_up_threshold", progress_deviative_up_threshold_, 1000.0);
    nh_.param<double>("progress_deviative_low_threshold", progress_deviative_low_threshold_, 500.0);
    nh_.param<double>("rc_deviative_up_threshold", rc_deviative_up_threshold_, 1000.0);
    nh_.param<double>("rc_deviative_low_threshold", rc_deviative_low_threshold_, 500.0);
    nh_.param<double>("init_yaw", init_yaw_, 0.0);
    nh_.param<double>("rc_expect_t", rc_expect_t_, 0.5);
    nh_.param<bool>("if_Auto_Yaw", if_Auto_Yaw_, true);
    nh_.param<int>("N_1", N_1_, 5.0);

    nh_.param<std::string>("hpoly_timeserial_path", hpoly_timeserial_path_, std::string("/offlines/HpolyTimeSerial.bin"));
    nh_.param<std::string>("glb_traj_path", glb_traj_path_, std::string("/offlines/global_trajectory.bin"));
    nh_.param<std::string>("corridor_path", corridor_path_, std::string("/offlines/corridors.bin"));

    std::cout << "[plan_manage] glb_traj_path_: " << glb_traj_path_ << std::endl;
    std::cout << "[plan_manage] corridor_path_: " << corridor_path_ << std::endl;
    std::cout << "[plan_manage] hpoly_timeserial_path_: " << hpoly_timeserial_path_ << std::endl;

    std::vector<double> tmp_gate_pos_vec;
    nh_.getParam("gate_pos_vec", tmp_gate_pos_vec);
    for (int i = 0; i < int(tmp_gate_pos_vec.size()) / 3; i++)
    {
      Eigen::Vector3d tmp_position;
      tmp_position << tmp_gate_pos_vec[3 * i], tmp_gate_pos_vec[3 * i + 1], tmp_gate_pos_vec[3 * i + 2];
      gate_pos_list_.push_back(tmp_position);
    }

    std::vector<double> tmp_gate_oritation;
    nh_.getParam("gate_oritation", tmp_gate_oritation);
    for (int i = 0; i < int(tmp_gate_oritation.size()); i++)
    {
      gate_oritation_vector.push_back(Eigen::Vector3d(cos(tmp_gate_oritation[i]), sin(tmp_gate_oritation[i]), 0));
    }

    terminal_pos_ = gate_pos_list_.back();
    terminal_vec_ = gate_oritation_vector.back();

    std::vector<double> initial_pos;
    nh_.getParam("initial_pos", initial_pos);
    init_x_ = initial_pos[0];
    init_y_ = initial_pos[1];
    init_z_ = initial_pos[2];

    std::vector<double> tmp_target_p;
    nh_.getParam("target_pos", tmp_target_p);
    global_target_p_ << tmp_target_p[0], tmp_target_p[1], tmp_target_p[2];
    global_yaw_target_p_ = global_target_p_ + 2.0 * terminal_vec_;
    lcYawTrajOptPtr_->setYawTarget(global_yaw_target_p_);

    // TODO: delete
    tmax = 0.0;

    traj_state_ = ASSIST;
    last_traj_state_ = traj_state_;

    current_id_ = 0;
    gate_index_ = 0;

    t_risk_trigger_ = ros::Time::now();

    // raed the msgs from binary
    ROS_WARN("[plan_manage] Start loding global traj and corridors !");
    ReadCorridors();
    SetEmergencyStopBox();
    ReadGlobalTraj();
    ROS_WARN("[PlanManager] initialization successfully !");
  }

  void PlanManager::SetEmergencyStopBox()
  {
    Eigen::MatrixXd safe_box = CorridorPtr_->getHpoly(CorridorPtr_->getCorridorNum() - 1);
    EmergencyStopOptPtr_->setSafeBox(safe_box);
  }

  void PlanManager::PrintState(Traj_STATE state)
  {
    std::string state_str[5] = {"PRIM", "ASSIST"};
    std::cout << "CURRENT FLY STATE:" << state_str[int(state)] << std::endl;
  }

  void PlanManager::StateChange(Traj_STATE state)
  {

    if (state != last_traj_state_)
    {
      std::string state_str[5] = {"PRIM", "ASSIST"};
      std::cout << "FLY STATE CHANGED: From " << state_str[int(last_traj_state_)] << " --> " << state_str[int(state)] << std::endl;
    }
    traj_state_ = state;
    last_traj_state_ = traj_state_;
  }

  void PlanManager::polyhedra_callback(const decomp_ros_msgs::PolyhedronArrayPtr &msgPtr)
  {
    if (receive_hpolys_)
      return;
    TunnelLikeHpolys_.clear();
    for (auto poly : msgPtr->polyhedrons)
    {
      Eigen::MatrixXd hp;
      int nump = poly.points.size();
      hp.resize(6, nump);
      for (int j = 0; j < nump; j++)
      {
        hp.col(j) << poly.normals[j].x, poly.normals[j].y, poly.normals[j].z, poly.points[j].x, poly.points[j].y, poly.points[j].z;
      }
      TunnelLikeHpolys_.push_back(hp);
    }
    CorridorPtr_->setCorridor(TunnelLikeHpolys_);
    receive_hpolys_ = true;
  }

  void PlanManager::GlobalTraj_callback(const quadrotor_msgs::PolynomialTrajectoryPtr &msgPtr)
  {
    if (get_glb_traj_)
      return;
    std::vector<CoefficientMat> coeflist;
    std::vector<double> times;

    int shift = 0;
    for (int idx = 0; idx < (int)msgPtr->num_segment; idx++)
    {
      int order = msgPtr->order[idx];
      CoefficientMat coefmat;
      coefmat.setZero();

      for (int j = 0; j <= order; j++)
      {
        coefmat(0, j) = msgPtr->coef_x[shift + j];
        coefmat(1, j) = msgPtr->coef_y[shift + j];
        coefmat(2, j) = msgPtr->coef_z[shift + j];
      }
      times.push_back(msgPtr->time[idx]);
      coeflist.push_back(coefmat);
      shift += (order + 1);
    }
    Trajectory traj(times, coeflist);
    glb_traj_ = traj;

    glb_start_timestamp_ = msgPtr->header.stamp;
    get_glb_traj_ = true;

    ROS_WARN("get global traj !");
    return;
  }

  void PlanManager::ReadGlobalTraj()
  {
    quadrotor_msgs::PolynomialTrajectory Trajmsg;
    wr_msg::readMsg(Trajmsg, ros::package::getPath("plan_manage") + glb_traj_path_);
    std::vector<CoefficientMat> coeflist;
    std::vector<double> times;

    int shift = 0;
    for (int idx = 0; idx < (int)Trajmsg.num_segment; idx++)
    {
      int order = Trajmsg.order[idx];
      CoefficientMat coefmat;
      coefmat.setZero();

      for (int j = 0; j <= order; j++)
      {
        coefmat(0, j) = Trajmsg.coef_x[shift + j];
        coefmat(1, j) = Trajmsg.coef_y[shift + j];
        coefmat(2, j) = Trajmsg.coef_z[shift + j];
      }
      times.push_back(Trajmsg.time[idx]);
      coeflist.push_back(coefmat);
      shift += (order + 1);
    }
    Trajectory traj(times, coeflist);
    glb_traj_ = traj;

    get_glb_traj_ = true;

    assistProgressOptPtr_->setRefTrajectory(glb_traj_);

    ROS_WARN("load global traj !");
  }

  void PlanManager::ReadCorridors()
  {
    decomp_ros_msgs::PolyhedronArray PolyhedronMsg;
    TunnelLikeHpolys_.clear();
    wr_msg::readMsg(PolyhedronMsg, ros::package::getPath("plan_manage") + corridor_path_);
    for (auto poly : PolyhedronMsg.polyhedrons)
    {
      Eigen::MatrixXd hp;
      int nump = poly.points.size();
      hp.resize(6, nump);
      for (int j = 0; j < nump; j++)
      {
        hp.col(j) << poly.normals[j].x, poly.normals[j].y, poly.normals[j].z, poly.points[j].x, poly.points[j].y, poly.points[j].z;
      }
      TunnelLikeHpolys_.push_back(hp);
    }

    CorridorPtr_->setCorridor(TunnelLikeHpolys_);

    quadrotor_msgs::HpolyTimeSerial timeSerial;
    wr_msg::readMsg(timeSerial, ros::package::getPath("plan_manage") + hpoly_timeserial_path_);

    std::cout << "TunnelLikeHpolys_.size()" << TunnelLikeHpolys_.size() << std::endl;
    std::cout << "timeSerial.t_start.size()" << timeSerial.t_start.size() << std::endl;

    for (int i = 0; i < (int)timeSerial.t_start.size(); i++)
    {
      time_serial_.push_back(std::make_pair(timeSerial.t_start[i], timeSerial.t_end[i]));
    }

    assistProgressOptPtr_->setTimeSerial(time_serial_);

    assert(timeSerial.t_start.size() == TunnelLikeHpolys_.size());

    receive_hpolys_ = true;

    ROS_WARN("load global corridors !");
  }

  void PlanManager::SafeTrajCallback(const ros::TimerEvent &event)
  {
    if (!circle_pos_receive_)
      return;
    int nearest_gate_index = 0;
    double min_dist = 999.0, dist_temp;

    for (int i = gate_index_; i < (int)gate_pos_list_.size(); i++)
    {
      dist_temp = (gate_pos_list_[i] - circle_pos_).norm();
      if (dist_temp < min_dist)
      {
        nearest_gate_index = i;
        min_dist = dist_temp;
      }
    }

    gate_pos_list_[nearest_gate_index] = circle_pos_;
    local_replan(nearest_gate_index);
  }

  void PlanManager::TicCallback(const ros::TimerEvent &event)
  {
    publish_desired_yaw();

    if ((odom_p_ - terminal_pos_).norm() < 1.0 && (odom_p_ - terminal_pos_).dot(terminal_vec_) > DBL_EPSILON && !reach_flag_ && gate_index_ > (int)gate_pos_list_.size() - 2)
    {
      std::cout << "[TicCallback] trigger_odom_p: " << odom_p_.transpose() << std::endl;
      static ros::Time end_timer = ros::Time::now();
      cout << "reach the target, the complete time: " << (end_timer - start_timer).toSec() << endl;
      Debug_msg_.finish_time = (end_timer - start_timer).toSec();
      cout << "tmax: " << tmax * 1000 << " ms " << endl;
      reach_flag_ = true;
    }

    static ros::Time t_vis_corridor = ros::Time::now();
    visualize_traj();
    if (receive_hpolys_ && (ros::Time::now() - t_vis_corridor).toSec() > 3.0)
    {
      t_vis_corridor = ros::Time::now();
      visCorridor(TunnelLikeHpolys_);
    }
  }

  void PlanManager::rcvOdometryCallback(const nav_msgs::Odometry &odom)
  {
    odom_p_(0) = odom.pose.pose.position.x;
    odom_p_(1) = odom.pose.pose.position.y;
    odom_p_(2) = odom.pose.pose.position.z;

    odom_v_(0) = odom.twist.twist.linear.x;
    odom_v_(1) = odom.twist.twist.linear.y;
    odom_v_(2) = odom.twist.twist.linear.z;
  }

  void PlanManager::rcvCircleCallback(const geometry_msgs::PoseStamped &msg)
  {
    circle_pos_receive_ = true;
    circle_pos_ << msg.pose.position.x, msg.pose.position.y, msg.pose.position.z;
  }

  void PlanManager::rcvJoyCmdCallback(const mavros_msgs::RCIn::ConstPtr &msg)
  {
    Eigen::VectorXi channels;
    int chansize = msg->channels.size();
    channels.resize(chansize);
    for (int i = 0; i < chansize; i++)
    {
      channels[i] = msg->channels[i];
    }

    // set trigger time
    RcPtr_->Inputs(channels);
    static bool trgger_flag = false;
    if (channels[5] == 1000 && !trgger_flag)
    {
      start_timer = ros::Time::now();
      glb_start_timestamp_ = ros::Time::now();
      trgger_flag = true;
      fly_start_trgger_ = true;
      ROS_WARN("flying start!");
    }

    if (!get_glb_traj_ || !fly_start_trgger_ || reach_flag_)
    {
      return;
    }

    // caculate the direct traj from the RC
    // for primitive only
    Trajectory7 prim_tmp_traj;
    Eigen::MatrixXd prim_initState(3, 4);
    ros::Time replan_timestamp;
    double replan_t;

    get_replan_state_prim(replan_timestamp, prim_initState, prim_tmp_traj, replan_t);
    Eigen::Vector3d RC_vel = RcPtr_->GetJoyVelCmd();
    RC_prim = calculateVelBasedPrimitive(RC_vel, prim_initState.col(0), prim_initState.col(1), prim_initState.col(2), prim_initState.col(3), rc_expect_t_);
    prim_tmp_traj.append(RC_prim);

    visPtr_->visualize_traj(prim_tmp_traj, "RC_Primitive");
    lc_prim_traj_ = prim_tmp_traj;
    lc_prim_start_timestamp_ = replan_timestamp;
    get_lc_prim_traj_ = true;
  }

  void PlanManager::local_replan(const int &gate_index)
  {
    if (!get_glb_traj_ || !get_lc_prim_traj_ || !receive_hpolys_ || !fly_start_trgger_)
      return;

    Debug_msg_.recover = 0.0;
    if (current_id_ == -1)
    {
      current_id_ = CorridorPtr_->getCurrentHpolyForward(odom_p_, 0);
      ROS_WARN("out of the hpolys, try to recover from current pos");
      Debug_msg_.recover = 1.0;
      Debug_pub.publish(Debug_msg_);
      return;
    }

    if (reach_flag_)
    {
      if (!if_emergency_stop_)
      {
        ros::Time emergency_t_start = ros::Time::now();
        ros::Time replan_timestamp_emergency;
        Trajectory tmp_traj_emergency;
        Eigen::MatrixXd initState_m(3, 3), finState_m(3, 3);
        double replan_t_emergency;

        get_replan_state(replan_timestamp_emergency, initState_m, finState_m, tmp_traj_emergency, gate_index, replan_t_emergency);

        Eigen::MatrixXd waypoints_emergency(3, 1);

        finState_m.col(0) << -6.0, 0, 1.0;
        finState_m.col(1).setZero();
        finState_m.col(2).setZero();

        waypoints_emergency = 0.5 * (initState_m.col(0) + finState_m.col(0));

        double dt_emergency = 1.0;

        Trajectory traj_emergency;
        double res_emergency = EmergencyStopOptPtr_->generate_traj(initState_m, finState_m, waypoints_emergency, dt_emergency, traj_emergency);

        if (res_emergency)
        {
          tmp_traj_emergency.append(traj_emergency);
          traj_emergency = tmp_traj_emergency;

          traj_emergency_ = traj_emergency;
          replan_timestamp_emergency_ = replan_timestamp_emergency;

          publish_traj(traj_emergency_, replan_timestamp_emergency);
          visPtr_->visualize_traj(traj_emergency_, "traj_emergency");

          if_emergency_stop_ = true;
          Debug_msg_.emergency_cal_time = (ros::Time::now() - emergency_t_start).toSec() * 1000;
          Debug_pub.publish(Debug_msg_);
        }
      }
      return;
    }
    //------------------------------------------------------PART1：get init state------------------------------------------------
    ros::Time t1, t2, t3, t4;
    t1 = ros::Time::now();

    ros::Time replan_timestamp;
    Trajectory7 tmp_traj;
    Trajectory1D tmp_yaw_traj;
    Eigen::MatrixXd initState(3, 4), finState(3, 4);
    Eigen::VectorXd yaw_initState(3), yaw_finState(3);
    double replan_t;
    double yaw_replan_t;
    get_replan_state(replan_timestamp, initState, finState, tmp_traj, gate_index, replan_t);
    int heuristic_current_id = current_id_ - 1 > 0 ? current_id_ - 1 : 0;
    current_id_ = CorridorPtr_->getCurrentHpolyForward(initState.col(0), heuristic_current_id);

    if (current_id_ == -1)
    {
      StateChange(ASSIST);
      return;
    }

    //------------------------------------------------------PART2：set the trajopt variables------------------------------------------------
    // TODO: if delete
    RiskPtr_->setCurrentState(initState);

    assistProgressOptPtr_->setCurrentCorridorId(current_id_);
    RiskPtr_->setCurrentCorridorId(current_id_);

    /*------------------------------------------------------PART3: get init waypoints and time for lc plan origin uniform use----------------------------------------*/
    int N_1 = N_1_;

    Eigen::MatrixXd waypoints(3, N_1 - 1); // for local traj
    Eigen::VectorXd Twaypoints(N_1 - 1);
    Trajectory7 cruise_traj;
    ros::Time cruise_traj_start_timestamp;

    double cruise_replan_time, t_terminal;
    double dt = replan_horizon_ / N_1;
    t_terminal = replan_horizon_;

    Debug_msg_.waypoint_check = 0.0;
    if (get_lc_traj_)
    {
      cruise_traj = lc_traj_;
      cruise_traj_start_timestamp = lc_start_timestamp_;
      cruise_replan_time = (replan_timestamp - cruise_traj_start_timestamp).toSec() + replan_duration_;

      for (int i = 0; i < N_1 - 1; i++)
      {
        double t_tmp = cruise_replan_time + (i + 1) * dt;
        if (t_tmp > t_terminal)
        {
          Eigen::MatrixXd terminal_state(3, 3);
          terminal_state << cruise_traj.getPos(t_terminal), cruise_traj.getVel(t_terminal), cruise_traj.getAcc(t_terminal);
          Eigen::MatrixXd state_tmp = StateIntegral(terminal_state, t_tmp - t_terminal);
          waypoints.col(i) = state_tmp.col(0);
          Debug_msg_.waypoint_check += 1.0;
        }
        else
        {
          waypoints.col(i) = cruise_traj.getPos(t_tmp);
        }
      }
    }
    else
    {
      cruise_replan_time = replan_duration_;
      for (int i = 0; i < N_1 - 1; i++)
      {
        waypoints.col(i) = initState.col(0);
      }
    }

    // set init Twaypoints
    Eigen::VectorXd initT(3), finT(3);
    Trajectory1D tmp_time_traj;
    YawBoundaryCond boundary;

    if (get_lc_traj_)
    {
      double tpos = lc_time_traj_.getPos(cruise_replan_time);
      double tvel = lc_time_traj_.getVel(cruise_replan_time);
      double tacc = lc_time_traj_.getAcc(cruise_replan_time);
      initT << tpos, tvel, tacc;
      boundary << lc_time_traj_.getPos(cruise_replan_time - replan_duration_),
          lc_time_traj_.getVel(cruise_replan_time - replan_duration_),
          lc_time_traj_.getAcc(cruise_replan_time - replan_duration_),
          tpos,
          tvel,
          tacc;
      tmp_time_traj = Trajectory1D(std::vector<Piece1D>{Piece1D(boundary, replan_duration_)});
    }
    else
    {
      initT << cruise_replan_time, 0, 0;
      boundary << cruise_replan_time - replan_duration_, 0, 0, cruise_replan_time, 0, 0;
      tmp_time_traj = Trajectory1D(std::vector<Piece1D>{Piece1D(boundary, replan_duration_)});
    }

    if (get_lc_traj_)
    {
      double t_end, tdot_end, tddot_end;
      for (int i = 0; i < N_1; i++)
      {
        double t_tmp = cruise_replan_time + (i + 1) * dt;
        if (i == N_1 - 1)
        {
          if (t_tmp > t_terminal)
          {
            t_end = lc_time_traj_.getPos(t_terminal);
            tdot_end = lc_time_traj_.getVel(t_terminal);
            tddot_end = lc_time_traj_.getAcc(t_terminal);
            finT << t_end + (t_tmp - t_terminal) * tdot_end + 0.5 * (t_tmp - t_terminal) * (t_tmp - t_terminal) + tddot_end,
                tdot_end + (t_tmp - t_terminal) * tddot_end,
                tddot_end;
          }
          else
          {
            finT << lc_time_traj_.getPos(t_tmp), lc_time_traj_.getVel(t_tmp), lc_time_traj_.getAcc(t_tmp);
          }
        }
        else
        {
          if (t_tmp > t_terminal)
          {
            t_end = lc_time_traj_.getPos(t_terminal);
            tdot_end = lc_time_traj_.getVel(t_terminal);
            tddot_end = lc_time_traj_.getAcc(t_terminal);
            Twaypoints(i) = t_end + (t_tmp - t_terminal) * tdot_end + 0.5 * (t_tmp - t_terminal) * (t_tmp - t_terminal) + tddot_end;
            Debug_msg_.time_check += 1.0;
          }
          else
          {
            Twaypoints(i) = lc_time_traj_.getPos(t_tmp);
          }
        }
      }
    }
    else
    {
      for (int i = 0; i < N_1 - 1; i++)
      {
        Twaypoints(i) = cruise_replan_time + (i + 1) * dt;
      }
      finT << cruise_replan_time + N_1 * dt, 0, 0;
    }

    t2 = ros::Time::now();
    assistProgressOptPtr_->setInitialState(initState, finState, initT, finT, N_1);
    assistProgressOptPtr_->setVelCmd(RcPtr_->GetJoyVelCmd());

    //------------------------------------------------------PART4：start local replan------------------------------------------------

    bool res = assistProgressOptPtr_->generate_traj(waypoints, Twaypoints, replan_horizon_, lc_traj_, lc_time_traj_, RC_prim);
    t3 = ros::Time::now();

    std::vector<double> Times = assistProgressOptPtr_->getComputationTimes();
    Debug_msg_.replan_traj_t_1 = (t2 - t1).toSec() * 1000;
    Debug_msg_.replan_traj_t_2 = (t3 - t2).toSec() * 1000;
    Debug_msg_.traj_generate_t_1 = Times[0];
    Debug_msg_.traj_generate_t_2 = Times[1];
    Debug_msg_.traj_generate_t_3 = Times[2];
    Debug_msg_.traj_generate_t_4 = Times[3];
    Debug_msg_.traj_generate_iter = assistProgressOptPtr_->getIteration();

    Debug_msg_.rhoTProgress = assistProgressOptPtr_->getRhoTProgress();

    double CorridorPnt = assistProgressOptPtr_->getCorridorPenalty();
    Debug_msg_.CorridorPnt = CorridorPnt;

    // debug used
    std::vector<Eigen::MatrixXd> HpolysIn = assistProgressOptPtr_->getHpolysIn();
    // visCorridor(HpolysIn);

    if (res)
    {
      tmp_traj.append(lc_traj_);
      lc_traj_ = tmp_traj;
      tmp_time_traj.append(lc_time_traj_);

      lc_time_traj_ = tmp_time_traj;
      lc_start_timestamp_ = replan_timestamp;
      gate_index_ = gate_index;
      get_lc_traj_ = true;

      publish_traj(lc_traj_, lc_start_timestamp_);

      visPtr_->visualize_1Dtraj(lc_time_traj_, "lc_time_traj");
    }

    //------------------------------------------------------PART5：set the local yaw plan------------------------------------------------
    if (!if_Auto_Yaw_)
    {
      std_msgs::Float32 t_opt;

      t_opt.data = (t3 - t1).toSec() * 1000;
      OptTime_pub_.publish(t_opt);

      // setDeviateCost
      if ((t3 - t1).toSec() > tmax)
      {
        tmax = (t3 - t1).toSec();
      }

      Debug_pub.publish(Debug_msg_);

      ROS_WARN("Replan sucessfully !");

      return;
    }

    get_yaw_replan_state(replan_timestamp, yaw_initState, yaw_finState, tmp_yaw_traj, gate_index, yaw_replan_t); // get yaw initial state

    t2 = ros::Time::now();
    bool yaw_res = false;

    // RC calc
    double RC_yaw_dot = RcPtr_->GetYawCmd();
    Yaw_prim = calculateYawPrimitive(yaw_initState, RC_yaw_dot, rc_expect_t_);

    // add a piece
    double end_p, end_v;
    double duration = Yaw_prim.getTotalDuration();
    end_p = Yaw_prim.getPos(duration);
    end_v = Yaw_prim.getVel(duration);
    YawCoefficientMat coeff;
    coeff << 0, 0, 0, 0, end_v, end_p;
    Trajectory1D aux_traj = Trajectory1D(std::vector<Piece1D>{Piece1D(3, coeff)});
    Yaw_prim.append(aux_traj);

    if (true)
    {

      double tmp_t = (replan_timestamp - lc_start_timestamp_).toSec() + replan_duration_;
      lcYawTrajOptPtr_->setForwardT(tmp_t);
      Eigen::Vector3d yaw_vec3D;
      Eigen::Vector2d yaw_vec;
      Eigen::VectorXd durations = lc_traj_.getDurations();
      Eigen::VectorXd D1 = durations.middleRows(1, N_1);
      Eigen::VectorXd yaw_waypoints(N_1); // for local traj
      Eigen::Vector3d vel_on_traj;
      double pieceD1 = (D1.sum() - (replan_timestamp - lc_start_timestamp_).toSec()) / N_1;

      if (pieceD1 < 0.01)
      {
        ROS_WARN("pieceD1 is too small! give up the yaw plan");
        return;
      }
      for (int i = 0; i < N_1; i++)
      {
        durations[i + 1] = pieceD1;
        D1[i] = pieceD1;
      }
      /*----------------- for local traj [ pos version] -------------------------*/
      Eigen::Vector3d pos_on_traj;
      for (int i = 0; i <= N_1; i++)
      {
        if (i < N_1)
        {
          pos_on_traj = lc_traj_.getPos(tmp_t);
          yaw_waypoints(i) = getYawVecToGate(pos_on_traj, gate_index);
          tmp_t += durations[i + 1];
        }
        else if (i == N_1 - 1)
        {
          pos_on_traj = lc_traj_.getPos(tmp_t);
          double angle_gate, angle_back;

          angle_gate = getYawVecToGate(pos_on_traj, gate_index);

          double tmp_tmp_t = tmp_t - 0.1;
          pos_on_traj = lc_traj_.getPos(tmp_tmp_t);

          angle_back = getYawVecToGate(pos_on_traj, gate_index);

          roundAngle(angle_gate, angle_back);

          yaw_finState << angle_gate, (angle_gate - angle_back) / 0.1, 0.01;
        }
      }

      roundAngle(yaw_waypoints[0], yaw_initState[0]);
      // reconstruct
      for (int i = 1; i < N_1; i++)
      {
        roundAngle(yaw_waypoints[i], yaw_waypoints[i - 1]);
      }
      roundAngle(yaw_finState[0], yaw_waypoints[N_1 - 1]);

      lcYawTrajOptPtr_->setOnTraj(lc_traj_);
      lcYawTrajOptPtr_->setAssistTraj(Yaw_prim);

      lcYawTrajOptPtr_->setGatesPosAndVec(gate_pos_list_, gate_oritation_vector);
      lcYawTrajOptPtr_->setGateIndex(gate_index);

      Eigen::VectorXd wp1 = yaw_waypoints.middleRows(1, N_1 - 1);

      t3 = ros::Time::now();
      yaw_res = lcYawTrajOptPtr_->generate_traj(yaw_initState, yaw_finState, wp1, D1, lc_yaw_traj_);
      t4 = ros::Time::now();

      Debug_msg_.replan_yaw_t_1 = (t3 - t2).toSec() * 1000;
      Debug_msg_.replan_yaw_t_2 = (t4 - t3).toSec() * 1000;
    }

    if (yaw_res)
    {
      tmp_yaw_traj.append(lc_yaw_traj_);
      visPtr_->visualize_1Dtraj(tmp_yaw_traj, "debug_yaw_traj");

      lc_yaw_traj_ = tmp_yaw_traj;

      lc_yaw_start_timestamp_ = replan_timestamp;
      get_lc_yaw_traj_ = true;
    }

    std_msgs::Float32 t_opt;

    t_opt.data = (t4 - t1).toSec() * 1000;
    OptTime_pub_.publish(t_opt);

    // setDeviateCost
    if ((t4 - t1).toSec() > tmax)
    {
      tmax = (t4 - t1).toSec();
    }

    Debug_pub.publish(Debug_msg_);

    ROS_WARN("Replan sucessfully !");
  }

  double PlanManager::getYawVecToGate(const Eigen::Vector3d pos, const int current_ids)
  {
    int id_tmp = -1;
    Eigen::Vector3d gate_pos;
    double angle;
    Eigen::Vector3d yaw_vec3D;
    Eigen::Vector2d yaw_vec;

    for (int j = current_ids; j < (int)gate_pos_list_.size(); j++)
    {
      gate_pos = gate_pos_list_[j];
      if ((gate_pos - pos).dot(gate_oritation_vector[j]) > 1e-2 && (gate_pos - pos).norm() > 1e-1)
      {
        yaw_vec3D = gate_pos - pos;
        yaw_vec = yaw_vec3D.head(2);
        angle = atan2(yaw_vec[1], yaw_vec[0] + 1e-6);
        id_tmp = 1;
        return angle;
      }
    }
    if (id_tmp == -1)
    {
      yaw_vec3D = global_yaw_target_p_ - pos;
      yaw_vec = yaw_vec3D.head(2);
      angle = atan2(yaw_vec[1], yaw_vec[0] + 1e-6);
      return angle;
    }
    ROS_ERROR("wrong angle calculated in init set !");
    return -100;
  }

  void PlanManager::get_replan_state(ros::Time &replan_timestamp, Eigen::MatrixXd &initState, Eigen::MatrixXd &finState, Trajectory &tmp_traj, const int &gate_index, double &replan_t)
  {
    replan_timestamp = ros::Time::now();
    Trajectory7 cruise_traj;
    ros::Time cruise_traj_start_timestamp;

    if (get_lc_traj_)
    {
      cruise_traj = lc_traj_;
      cruise_traj_start_timestamp = lc_start_timestamp_;
      replan_t = (replan_timestamp - cruise_traj_start_timestamp).toSec() + replan_duration_;
      double replan_terminal = replan_t + replan_horizon_;
      double cruise_terminal = replan_horizon_;
      double t_tmp = replan_terminal < cruise_terminal ? replan_terminal : cruise_terminal;
      target_p_ = cruise_traj.getPos(t_tmp);
      target_v_ = cruise_traj.getVel(t_tmp);
      target_a_ = cruise_traj.getAcc(t_tmp);
      target_j_ = cruise_traj.getJerk(t_tmp);

      if (replan_terminal > cruise_terminal + 1e-6)
      {
        finState.col(2) = target_a_ + replan_t * target_j_;
        finState.col(1) = target_v_ + replan_t * target_a_ + 0.5 * replan_t * replan_t * target_j_;
        finState.col(0) = target_p_ + replan_t * target_v_ + 0.5 * replan_t * replan_t * target_a_ + 0.167 * replan_t * replan_t * replan_t * target_j_;
      }

      start_p_ = cruise_traj.getPos(replan_t);
      start_v_ = cruise_traj.getVel(replan_t);
      start_a_ = cruise_traj.getAcc(replan_t);

      initState << start_p_, start_v_, start_a_;

      BoundaryCond boundary;
      boundary << cruise_traj.getPos(replan_t - replan_duration_),
          cruise_traj.getVel(replan_t - replan_duration_),
          cruise_traj.getAcc(replan_t - replan_duration_),
          cruise_traj.getPos(replan_t),
          cruise_traj.getVel(replan_t),
          cruise_traj.getAcc(replan_t);

      tmp_traj = Trajectory(std::vector<Piece>{Piece(boundary, replan_duration_)});
    }
    else
    {
      replan_t = replan_timestamp.toSec() + replan_duration_;
      start_p_ << init_x_, init_y_, init_z_;
      start_v_ << 0, 0, 0;
      start_a_ << 0, 0, 0;

      initState << start_p_, start_v_, start_a_;

      finState = initState;

      BoundaryCond boundary;
      boundary << start_p_, start_v_, start_a_, start_p_, start_v_, start_a_;

      tmp_traj = Trajectory(std::vector<Piece>{Piece(boundary, replan_duration_)});
    }

    return;
  }

  void PlanManager::get_replan_state(ros::Time &replan_timestamp, Eigen::MatrixXd &initState, Eigen::MatrixXd &finState, Trajectory7 &tmp_traj, const int &gate_index, double &replan_t)
  {
    replan_timestamp = ros::Time::now();
    Trajectory7 cruise_traj;
    ros::Time cruise_traj_start_timestamp;

    if (get_lc_traj_)
    {
      cruise_traj = lc_traj_;
      cruise_traj_start_timestamp = lc_start_timestamp_;
      replan_t = (replan_timestamp - cruise_traj_start_timestamp).toSec() + replan_duration_;
      double replan_terminal = replan_t + replan_horizon_;
      double cruise_terminal = replan_horizon_;
      double t_tmp = replan_terminal < cruise_terminal ? replan_terminal : cruise_terminal;

      target_p_ = cruise_traj.getPos(t_tmp);
      target_v_ = cruise_traj.getVel(t_tmp);
      target_a_ = cruise_traj.getAcc(t_tmp);
      target_j_ = cruise_traj.getJerk(t_tmp);

      if (replan_terminal > cruise_terminal + 1e-6)
      {
        finState.col(3) = target_j_;
        finState.col(2) = target_a_ + replan_t * target_j_;
        finState.col(1) = target_v_ + replan_t * target_a_ + 0.5 * replan_t * replan_t * target_j_;
        finState.col(0) = target_p_ + replan_t * target_v_ + 0.5 * replan_t * replan_t * target_a_ + 0.167 * replan_t * replan_t * replan_t * target_j_;
      }

      if (traj_state_ == PRIM)
      {
        cruise_traj = lc_prim_traj_;
        cruise_traj_start_timestamp = lc_prim_start_timestamp_;
        replan_t = (replan_timestamp - cruise_traj_start_timestamp).toSec() + replan_duration_;
      }

      start_p_ = cruise_traj.getPos(replan_t);
      start_v_ = cruise_traj.getVel(replan_t);
      start_a_ = cruise_traj.getAcc(replan_t);
      start_j_ = cruise_traj.getJerk(replan_t);
      initState << start_p_, start_v_, start_a_, start_j_;

      BoundaryCond7 boundary;
      boundary << cruise_traj.getPos(replan_t - replan_duration_),
          cruise_traj.getVel(replan_t - replan_duration_),
          cruise_traj.getAcc(replan_t - replan_duration_),
          cruise_traj.getJerk(replan_t - replan_duration_),
          start_p_,
          start_v_,
          start_a_,
          start_j_;

      tmp_traj = Trajectory7(std::vector<Piece7>{Piece7(boundary, replan_duration_)});
    }
    else
    {
      replan_t = replan_timestamp.toSec() + replan_duration_;
      start_p_ = odom_p_;
      start_v_ << 0, 0, 0;
      start_a_ << 0, 0, 0;
      start_j_ << 0, 0, 0;

      initState << start_p_, start_v_, start_a_, start_j_;
      finState = initState;

      BoundaryCond7 boundary;
      boundary << start_p_, start_v_, start_a_, start_j_,
          start_p_, start_v_, start_a_, start_j_;
      tmp_traj = Trajectory7(std::vector<Piece7>{Piece7(boundary, replan_duration_)});
    }

    return;
  }

  Eigen::MatrixXd PlanManager::StateIntegral(const Eigen::MatrixXd State, const double integralT)
  {
    Eigen::MatrixXd finState(3, 3);
    finState.col(2) = State.col(2);
    finState.col(1) = State.col(1) + integralT * State.col(2);
    finState.col(0) = State.col(0) + integralT * State.col(1) + 0.5 * integralT * integralT * State.col(2);
    return finState;
  }

  void PlanManager::get_yaw_replan_state(ros::Time &replan_timestamp, Eigen::VectorXd &initState, Eigen::VectorXd &finState, Trajectory1D &tmp_traj, const int &gate_index, double &replan_t)
  {
    Trajectory1D cruise_yaw_traj;
    ros::Time cruise_traj_start_timestamp;
    if (get_lc_yaw_traj_)
    {
      cruise_yaw_traj = lc_yaw_traj_;
      cruise_traj_start_timestamp = lc_yaw_start_timestamp_;

      replan_timestamp = ros::Time::now();
      replan_t = (replan_timestamp - cruise_traj_start_timestamp).toSec() + replan_duration_;

      initState << cruise_yaw_traj.getPos(replan_t), cruise_yaw_traj.getVel(replan_t), cruise_yaw_traj.getAcc(replan_t);

      YawBoundaryCond boundary;
      boundary << cruise_yaw_traj.getPos(replan_t - replan_duration_),
          cruise_yaw_traj.getVel(replan_t - replan_duration_),
          cruise_yaw_traj.getAcc(replan_t - replan_duration_),
          cruise_yaw_traj.getPos(replan_t),
          cruise_yaw_traj.getVel(replan_t),
          cruise_yaw_traj.getAcc(replan_t);

      tmp_traj = Trajectory1D(std::vector<Piece1D>{Piece1D(boundary, replan_duration_)});
    }
    else
    {
      replan_timestamp = ros::Time::now();
      replan_t = replan_duration_;

      initState << init_yaw_, 0, 0;

      YawBoundaryCond boundary;
      boundary << init_yaw_,
          0,
          0,
          init_yaw_,
          0,
          0;

      tmp_traj = Trajectory1D(std::vector<Piece1D>{Piece1D(boundary, replan_duration_)});
    }

    return;
  }

  void PlanManager::get_replan_state_prim(ros::Time &replan_timestamp, Eigen::MatrixXd &initState, Trajectory7 &tmp_traj, double &replan_t)
  {
    replan_timestamp = ros::Time::now();
    Trajectory7 cruise_traj;
    ros::Time cruise_traj_start_timestamp;
    if (get_lc_prim_traj_ && get_lc_traj_)
    {
      if (traj_state_ == PRIM)
      {
        cruise_traj = lc_prim_traj_;
        cruise_traj_start_timestamp = get_lc_prim_traj_ ? lc_prim_start_timestamp_ : glb_start_timestamp_;
      }
      else if (traj_state_ == ASSIST)
      {
        cruise_traj = lc_traj_;
        cruise_traj_start_timestamp = get_lc_traj_ ? lc_start_timestamp_ : glb_start_timestamp_;
      }

      replan_t = (replan_timestamp - cruise_traj_start_timestamp).toSec() + replan_duration_;

      start_p_ = cruise_traj.getPos(replan_t);
      start_v_ = cruise_traj.getVel(replan_t);
      start_a_ = cruise_traj.getAcc(replan_t);
      start_j_ = cruise_traj.getJerk(replan_t);

      initState << start_p_, start_v_, start_a_, start_j_;
      BoundaryCond7 boundary;
      boundary << cruise_traj.getPos(replan_t - replan_duration_),
          cruise_traj.getVel(replan_t - replan_duration_),
          cruise_traj.getAcc(replan_t - replan_duration_),
          cruise_traj.getJerk(replan_t - replan_duration_),
          start_p_,
          start_v_,
          start_a_,
          start_j_;

      tmp_traj = Trajectory7(std::vector<Piece7>{Piece7(boundary, replan_duration_)});
    }
    else
    {
      replan_t = replan_timestamp.toSec() + replan_duration_;
      start_p_ = odom_p_;
      start_v_ << 0, 0, 0;
      start_a_ << 0, 0, 0;
      start_j_ << 0, 0, 0;

      initState << start_p_, start_v_, start_a_, start_j_;

      BoundaryCond7 boundary;
      boundary << start_p_,
          start_v_,
          start_a_,
          start_j_,
          start_p_,
          start_v_,
          start_a_,
          start_j_;

      tmp_traj = Trajectory7(std::vector<Piece7>{Piece7(boundary, replan_duration_)});
    }

    return;
  }

  void PlanManager::visualize_traj()
  {
    if (get_glb_traj_)
    {
      visPtr_->visualize_traj(glb_traj_, "global_traj");
    }
    if (get_lc_traj_)
    {
      visPtr_->visualize_traj(lc_traj_, "local_traj");
    }
    if (get_lc_yaw_traj_)
    {
      visPtr_->visualize_1Dtraj(lc_yaw_traj_, "local_yaw_traj");
    }
  }

  void PlanManager::publish_traj(const Trajectory &traj, const ros::Time &traj_timestamp)
  {
    quadrotor_msgs::PolynomialTrajectory traj_msg;
    traj_msg = traj2msg(traj, traj_timestamp);
    traj_pub_.publish(traj_msg);
  }

  void PlanManager::publish_traj(const Trajectory7 &traj, const ros::Time &traj_timestamp)
  {
    quadrotor_msgs::PolynomialTrajectory traj_msg;
    traj_msg = traj2msg(traj, traj_timestamp);
    traj_pub_.publish(traj_msg);
  }

  quadrotor_msgs::PolynomialTrajectory PlanManager::traj2msg(const Trajectory &traj, const ros::Time &traj_timestamp)
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

  quadrotor_msgs::PolynomialTrajectory PlanManager::traj2msg(const Trajectory7 &traj, const ros::Time &traj_timestamp)
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
        CoefficientMat7 coemat = traj[i].getCoeffMat();
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

  void PlanManager::publish_desired_yaw()
  {
    if (!get_glb_traj_)
      return;
    std_msgs::Float32 desired_yaw;

    if (if_Auto_Yaw_)
    {
      ros::Time tnow = ros::Time::now();
      Trajectory1D cruise_yaw_traj;
      ros::Time cruise_start_time;
      double t;
      double des_yaw = init_yaw_;

      if (get_lc_yaw_traj_)
      {
        cruise_yaw_traj = lc_yaw_traj_;
        cruise_start_time = lc_yaw_start_timestamp_;

        t = (tnow - cruise_start_time).toSec();

        if (t >= cruise_yaw_traj.getTotalDuration())
          return;
        des_yaw = cruise_yaw_traj.getPos(t);

        // DEBUG
        Trajectory7 cruise_traj = lc_traj_;
        ros::Time cruise_traj_start_timestamp = lc_start_timestamp_;
        t = (tnow - cruise_traj_start_timestamp).toSec();
        Eigen::Vector3d pos = cruise_traj.getPos(t);
        Eigen::Vector3d vec(cos(des_yaw), sin(des_yaw), 0);
        visPtr_->visualize_arrow(pos, pos + vec, "yaw_vec");
      }

      desired_yaw.data = des_yaw;
      RcPtr_->SetYaw(des_yaw);
    }
    else
    {
      desired_yaw.data = RcPtr_->GetYaw();
    }

    desired_yaw_pub_.publish(desired_yaw);
  }

  Trajectory PlanManager::calculatePrimitive(Eigen::Vector3d av_cmd,
                                             Eigen::Vector3d q_pos,
                                             Eigen::Vector3d q_vel,
                                             double dur,
                                             bool isAcc)
  {
    Trajectory Primitive;
    CoefficientMat coeff;
    if (isAcc)
    {
      // add gravity
      Eigen::Vector3d acc_prim, gravity(0, 0, -9.8);
      Eigen::Vector3d rotor_drag;
      Eigen::Vector3d drag_coeff;
      drag_coeff << 0.4, 0.4, 0.2;
      rotor_drag = drag_coeff.asDiagonal() * q_vel;
      acc_prim = av_cmd + gravity - rotor_drag;
      coeff << Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 0.5 * acc_prim, q_vel, q_pos;
      Primitive = Trajectory(std::vector<Piece>{Piece(dur, coeff)});
    }
    else
    {
      coeff << Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), av_cmd, q_pos;
      Primitive = Trajectory(std::vector<Piece>{Piece(dur, coeff)});
    }

    return Primitive;
  }

  Trajectory7 PlanManager::calculateVelBasedPrimitive(Eigen::Vector3d v_cmd,
                                                      Eigen::Vector3d q_pos,
                                                      Eigen::Vector3d q_vel,
                                                      Eigen::Vector3d q_acc,
                                                      Eigen::Vector3d q_jrk,
                                                      double dur)
  {
    double a_max = config_.max_acc;
    double t_decelerate = (q_vel - v_cmd).norm() / a_max;
    t_decelerate = t_decelerate < 0.01 ? 0.01 : t_decelerate;
    Trajectory7 Primitive;
    int cnt = 0;
    double max_acc = 1000;
    while (max_acc > config_.max_acc)
    {
      if (cnt >= 40)
        break;
      BoundaryCond7 boundary;
      boundary << q_pos,
          q_vel,
          q_acc,
          q_jrk,
          v_cmd,
          Eigen::Vector3d::Zero(),
          Eigen::Vector3d::Zero(),
          Eigen::Vector3d::Zero();

      Primitive = Trajectory7(std::vector<Piece7>{Piece7(boundary, t_decelerate, true)});
      max_acc = Primitive.getMaxAccRate();
      cnt++;
      t_decelerate += 0.05;
    }
    return Primitive;
  }

  Trajectory PlanManager::calculateStopPrimitive(Eigen::Vector3d end_pos,
                                                 Eigen::Vector3d q_pos,
                                                 Eigen::Vector3d q_vel,
                                                 Eigen::Vector3d q_acc,
                                                 double dur)
  {
    Trajectory Primitive;
    CoefficientMat coeff;

    BoundaryCond boundary;
    boundary << q_pos,
        q_vel,
        q_acc,
        end_pos,
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Zero();

    Primitive = Trajectory(std::vector<Piece>{Piece(boundary, dur)});

    return Primitive;
  }

  Trajectory1D PlanManager::calculateYawPrimitive(const Eigen::Vector3d yaw_init,
                                                  const double yaw_dot_cmd,
                                                  const double dur)
  {
    Trajectory1D Primitive;
    YawCoefficientMat coeff;

    double yaw_end = (yaw_dot_cmd + yaw_init[1]) / 2 * dur + yaw_init[0];
    YawBoundaryCond boundary;
    boundary << yaw_init[0], yaw_init[1], yaw_init[2], yaw_end, yaw_dot_cmd, 0;

    Primitive = Trajectory1D(std::vector<Piece1D>{Piece1D(boundary, dur)});

    return Primitive;
  }

} // namespace planner