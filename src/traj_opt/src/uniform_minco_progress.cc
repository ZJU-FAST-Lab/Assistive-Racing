#include <traj_opt/uniform_minco_progress.h>

#include <traj_opt/chalk.hpp>
#include <optimization_opt/lbfgs_raw.hpp>

namespace traj_opt
{
  UniformMincoProgress::UniformMincoProgress(ros::NodeHandle &nh, Config &conf)
      : nh_(nh), config_(conf)
  {
    squared_vel_limit_ = config_.max_vel * config_.max_vel;
    squared_acc_limit_ = config_.max_acc * config_.max_acc;
    squared_jrk_limit_ = config_.max_jrk * config_.max_jrk;
    dist0_ = config_.dist0;
  }

  bool UniformMincoProgress::generate_traj(const Eigen::MatrixXd &Poswaypoints,
                                           const Eigen::VectorXd &Twaypoints,
                                           const double &Total_T,
                                           Trajectory7 &traj,
                                           Trajectory1D &time_traj,
                                           Trajectory7 &primitive)
  {
    ros::Time t1 = ros::Time::now();
    Prim = primitive;

    Eigen::Vector3d end_p, end_v;
    double duration_ = primitive.getTotalDuration();
    end_p = primitive.getPos(duration_);
    end_v = primitive.getVel(duration_);
    CoefficientMat7 coeff;
    coeff << Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), end_v, end_p;
    Trajectory7 aux_traj = Trajectory7(std::vector<Piece7>{Piece7(2, coeff)});
    Prim.append(aux_traj);

    waypoints_1_ = Poswaypoints;
    Twaypoints_1_ = Twaypoints;

    waypoints_dim_ = waypoints_1_.rows();
    dim_p_1_ = waypoints_dim_ * waypoints_1_.cols();

    dim_t_p_ = Twaypoints_1_.rows();

    state_order_ = end_state_.cols();
    dim_end_state_ = waypoints_dim_ * state_order_;
    dim_tend_state_ = fin_T_.size();

    total_t_ = Total_T;

    x_.resize(dim_p_1_ + dim_t_p_ + dim_end_state_ + dim_tend_state_);
    Eigen::Map<Eigen::MatrixXd> P1(x_.data(), waypoints_dim_, N_1_ - 1);
    Eigen::Map<Eigen::VectorXd> TP1(x_.data() + dim_p_1_, dim_t_p_);
    Eigen::Map<Eigen::MatrixXd> end_state(x_.data() + dim_p_1_ + dim_t_p_, waypoints_dim_, state_order_);
    Eigen::Map<Eigen::VectorXd> t_end_state(x_.data() + dim_p_1_ + dim_t_p_ + dim_end_state_, dim_tend_state_);

    P1 = waypoints_1_;
    TP1 = Twaypoints_1_;
    end_state = end_state_;
    t_end_state = fin_T_;

    t_lower_bound = config_.rc_t_lower_bound;
    K_fix_ = config_.rc_K_fix;
    step_fix_ = t_lower_bound / K_fix_;

    snapOpt_1_.generate(P1, total_t_);
    Timeopt_1_.generate(TP1, total_t_);

    ros::Time t2 = ros::Time::now();
    perperation_t1_ = (t2 - t1).toSec() * 1000;

    calculateSafeHpolys(snapOpt_1_.getTraj(), Timeopt_1_.getTraj());
    ros::Time t3 = ros::Time::now();
    perperation_t2_ = (t3 - t2).toSec() * 1000;

    calculatePrimSafeHpolys();
    ros::Time t4 = ros::Time::now();
    perperation_t3_ = (t4 - t3).toSec() * 1000;

    if (config_.dynamicRhoT)
      rhoTprogress_ = config_.rhoTProgress * vel_cmd_.norm() / config_.max_vel + 1;
    else
      rhoTprogress_ = config_.rhoTProgress;

    int opt_ret = optimize();

    ros::Time t5 = ros::Time::now();
    optimization_t_ = (t5 - t4).toSec() * 1000;

    if (opt_ret < 0)
    {
      return false;
    }

    snapOpt_1_.resetTailCon(end_state);
    snapOpt_1_.generate(P1, total_t_);
    Timeopt_1_.resetTailCon(t_end_state);
    Timeopt_1_.generate(TP1, total_t_);

    traj = snapOpt_1_.getTraj();
    time_traj = Timeopt_1_.getTraj();

    return true;
  }

  int UniformMincoProgress::optimize(const double &delta)
  {
    // Setup for L-BFGS solver
    lbfgs::lbfgs_parameter_t lbfgs_params;
    lbfgs_params.mem_size = 16;
    lbfgs_params.past = 3;
    lbfgs_params.g_epsilon = 0.0;
    lbfgs_params.min_step = 1e-32;
    lbfgs_params.delta = 1.0e-4;
    double minObjective;
    iteration_ = 0;
    auto ret = lbfgs::lbfgs_optimize(
        x_, minObjective, UniformMincoProgress::objectiveFunc, nullptr,
        UniformMincoProgress::earlyExit, this, lbfgs_params);
    return ret;
  }

  inline int UniformMincoProgress::earlyExit(void *ptrObj, const Eigen::VectorXd &x,
                                             const Eigen::VectorXd &grad, const double fx,
                                             const double step, int k, int ls)
  {
    // return k > 1e3;
    UniformMincoProgress *obj = reinterpret_cast<UniformMincoProgress *>(ptrObj);
    if (obj->iteration_ > 800)
      return 1;
    return 0;
  }

  void UniformMincoProgress::calculatePrimSafeHpolys()
  {
    d_set_.clear();
    double s = 0;
    Eigen::Vector3d pos_des;

    for (int i = 0; i <= K_fix_; i++)
    {
      pos_des = Prim.getPos(s);
      double d = CorridorPtr_->getBestDistInHpoly(pos_des, current_id_, config_.search_forward_idx);
      d_set_.push_back(d);
      s += step_fix_;
    }
  }

  void UniformMincoProgress::calculateSafeHpolys(const Trajectory7 pos_traj, const Trajectory1D time_traj)
  {
    Hpolys_in_.clear();
    std::vector<int> num_counter(CorridorPtr_->getCorridorNum());
    std::fill(num_counter.begin(), num_counter.end(), 0);

    int total_num_ball = N_1_ * config_.K + 1;
    int num_check = (N_1_ - 1) * config_.K;
    double dt = total_t_ / total_num_ball;
    double s = 0;
    Eigen::Vector3d pos, pos_ref, pos_on_glb;
    double t_ref;
    int id_tmp;
    int id_last = current_id_;
    double t_best = -1.0;

    ids_.clear();
    for (int i = 0; i < num_check; i++)
    {
      id_tmp = -1;
      pos = pos_traj.getPos(s);
      t_ref = time_traj.getPos(s);
      pos_on_glb = ref_traj_.getPos(t_ref);
      int heuristic_id = current_id_ - 1 >= 0 ? current_id_ - 1 : 0;
      if (i == 0)
      {
        pos_ref = findNeareastPos(pos, t_ref, current_id_, t_best);
      }
      else
      {
        pos_ref = findNeareastPos(pos, t_ref, id_last, t_best);
      }

      if (t_best == -1)
        t_best = t_ref;
      for (int k = heuristic_id; k < CorridorPtr_->getCorridorNum(); k++)
      {
        if (t_best >= time_searial_[k].first && t_best < time_searial_[k].second)
        {
          id_tmp = k;
          break;
        }
      }
      if (id_tmp == -1)
      {
        if (i == 0)
        {
          id_tmp = current_id_;
        }
        else
          id_tmp = ids_[i - 1];
      }

      if (id_tmp < id_last)
        id_tmp = id_last; // maybe important

      id_last = id_tmp;

      num_counter[id_tmp]++;
      ids_.push_back(id_tmp);
      s += dt;
    }

    for (int i = 0; i < CorridorPtr_->getCorridorNum(); i++)
    {
      if (num_counter[i] > 0)
      {
        Hpolys_in_.push_back(CorridorPtr_->getHpoly(i));
      }
    }
  }

  Eigen::Vector3d UniformMincoProgress::findNeareastPos(const Eigen::Vector3d pos, const double heuristic_t, int id_last, double &t_best)
  {
    double t_len = total_t_ * 0.5;
    double t_start = heuristic_t - t_len > 0 ? heuristic_t - t_len : 0;
    double t_end = heuristic_t + t_len < ref_traj_.getTotalDuration() ? heuristic_t + t_len : ref_traj_.getTotalDuration();
    double dist = 1000, dist_tmp, t_b = -1;
    // double t_best = 0;
    std::vector<Eigen::Vector3d> pos_set;
    std::vector<double> time_set;
    Eigen::Vector3d pos_best;
    Eigen::VectorXd durations = ref_traj_.getDurations();

    int idx_start = ref_traj_.locatePieceIdx(t_start);
    int idx_end = ref_traj_.locatePieceIdx(t_end);
    double t_accumulate = durations.head(idx_start).sum();

    for (int i = idx_start; i <= idx_end; i++)
    {
      Piece piece = ref_traj_.getPiece(i);
      if (i == idx_start)
        pos_set = piece.get_all_project_pt(pos, 0.0, time_set);
      else
        pos_set = piece.get_all_project_pt(pos, 0.0, time_set);
      if (!pos_set.empty())
      {
        for (int k = 0; k < int(pos_set.size()); k++)
        {
          auto pos_candidate = pos_set[k];
          if (CorridorPtr_->SafeCheckHpolyRange(pos_candidate, id_last, 5, 5) == -1)
            continue;
          dist_tmp = (pos - pos_candidate).squaredNorm();
          if (dist_tmp < dist)
          {
            dist = dist_tmp;
            pos_best = pos_candidate;
            t_b = time_set[k] + t_accumulate;
          }
        }
      }
      t_accumulate += piece.getDuration();
    }
    t_best = t_b;
    return pos_best;
  }

  Eigen::Vector3d UniformMincoProgress::findNeareastPos_numerical(const Eigen::Vector3d pos, const double heuristic_t, int id_last, double &t_best)
  {
    double dt = 0.01;
    double t_len = total_t_ * 0.5;
    double t_start = heuristic_t - t_len > 0 ? heuristic_t - t_len : 0;
    double t_end = heuristic_t + t_len < ref_traj_.getTotalDuration() ? heuristic_t + t_len : ref_traj_.getTotalDuration();
    double dist = 1000, dist_tmp;
    t_best = 0;

    for (double t = t_start; t < t_end; t += dt)
    {
      if (CorridorPtr_->SafeCheckHpolyRange(ref_traj_.getPos(t), id_last, 5, 5) == -1)
        continue;
      dist_tmp = (ref_traj_.getPos(t) - pos).squaredNorm();
      if (dist_tmp < 0)
      {
        std::cout << "wrong value: t: " << t << " ref_pos: " << ref_traj_.getPos(t) << " pos: " << pos << std::endl;
      }

      if (dist_tmp < dist)
      {
        dist = dist_tmp;
        t_best = t;
      }
    }
    return ref_traj_.getPos(t_best);
  }

  void UniformMincoProgress::setInitialState(const Eigen::MatrixXd &iniState,
                                             const Eigen::MatrixXd &finState,
                                             const Eigen::VectorXd &initT,
                                             const Eigen::VectorXd &finT,
                                             const int N)
  {
    ini_state_ = iniState;
    end_state_ = finState;
    ini_T_ = initT;
    fin_T_ = finT;

    N_1_ = N;

    Timeopt_2_.reset(N_1_);

    roundingState(ini_state_, ini_T_);
    roundingState(end_state_, fin_T_);
    setBoundConds(ini_state_, end_state_, ini_T_, fin_T_);
  }

  inline void UniformMincoProgress::roundingState(Eigen::MatrixXd &state, Eigen::VectorXd &T)
  {
    double tempNorm = state.col(1).norm();
    if (tempNorm > config_.max_vel)
    {
      state.col(1) *= (config_.max_vel - 0.001) / tempNorm;
      state.col(2) *= 0.0;
    }

    tempNorm = state.col(2).norm();
    if (tempNorm > config_.max_acc)
    {
      state.col(2) *= (config_.max_acc - 0.001) / tempNorm;
    }

    tempNorm = abs(T[1]);
    if (tempNorm > config_.T_max_vel)
    {
      T[1] *= (config_.T_max_vel - 0.001) / tempNorm;
    }

    tempNorm = abs(T[2]);
    if (tempNorm > config_.T_max_acc)
    {
      T[2] *= (config_.T_max_acc - 0.001) / tempNorm;
    }
  }

  inline void UniformMincoProgress::setBoundConds(const Eigen::MatrixXd &iniState,
                                                  const Eigen::MatrixXd &finState,
                                                  const Eigen::VectorXd &iniT,
                                                  const Eigen::VectorXd &finT)
  {
    // jerkOpt_1_.reset(iniState, finState, N_1_);
    snapOpt_1_.reset(iniState, finState, N_1_);
    Timeopt_1_.reset(iniT, finT, N_1_);
  }

  double UniformMincoProgress::calGatePenalty(const Eigen::MatrixXd &state, Eigen::Ref<Eigen::MatrixXd> grad)
  {
    double cost = 0.0;
    double violation = (state.col(0) - gate_pos_).squaredNorm();
    cost = config_.pnlGate * violation;
    grad.col(0) += config_.pnlGate * 2.0 * (state.col(0) - gate_pos_);
    return cost;
  }

  double UniformMincoProgress::calTimeIntPenalty(const int &N,
                                                 Eigen::MatrixX3d &gdC,
                                                 Eigen::VectorXd &gdCt)
  {
    Eigen::Vector3d pos, vel, acc, jrk;
    double t, t_dot, t_ddot, t_dddot;
    Eigen::Vector3d pos_ref, vel_ref;

    Eigen::Vector3d grad_tmp, grad_ref_tmp;
    double cost_tmp;
    double tgrad_tmp;
    double tcost_tmp;

    int idx = 0;

    static Eigen::Matrix<double, 6, 1> beta0, beta1, beta2, beta3;
    static Eigen::Matrix<double, 8, 1> snap_beta0, snap_beta1, snap_beta2, snap_beta3, snap_beta4;
    double s1, s2, s3, s4, s5, s6, s7;
    double step_vaj, omega;
    Eigen::Matrix<double, 8, 3> grad_Jv_to_c, grad_Ja_to_c, grad_Jj_to_c, grad_JCone_to_c;
    Eigen::Matrix<double, 8, 3> integral_grad_Jv_to_c, integral_grad_Ja_to_c, integral_grad_Jj_to_c, integral_grad_JCone_to_c;
    Eigen::Matrix<double, 6, 1> integral_grad_Jtv_to_ct, integral_grad_Jta_to_ct, integral_grad_Jtj_to_ct;

    Eigen::Matrix<double, 8, 3> integral_grad_Jprogress_to_c;
    Eigen::Matrix<double, 6, 1> integral_grad_Jprogress_to_ct, integral_grad_Jtprogress_to_ct, integral_grad_Jtpositive_to_ct, integral_grad_Jtbound_to_ct;

    vel_pnt_ = 0.0;
    acc_pnt_ = 0.0;
    jrk_pnt_ = 0.0;
    Cone_pnt_ = 0.0;
    t_progress_cost_ = 0.0;
    pos_progress_cost_ = 0.0;
    t_positive_pnt_ = 0.0;
    t_bound_pnt_ = 0.0;
    t_vel_pnt_ = 0.0;
    t_acc_pnt_ = 0.0;
    t_jrk_pnt_ = 0.0;

    grad_JCone_to_c.setZero();

    double ref_traj_duration = ref_traj_.getTotalDuration();
    double rho;

    for (int i = 0; i < N; ++i)
    {
      const auto &c = snapOpt_1_.getCoeffs().block<8, 3>(i * 8, 0);
      const auto &ct = Timeopt_1_.getCoeffs().block<6, 1>(i * 6, 0);
      integral_grad_Jv_to_c.setZero();
      integral_grad_Ja_to_c.setZero();
      integral_grad_Jj_to_c.setZero();
      integral_grad_JCone_to_c.setZero();
      integral_grad_Jprogress_to_c.setZero();

      integral_grad_Jtv_to_ct.setZero();
      integral_grad_Jta_to_ct.setZero();
      integral_grad_Jtj_to_ct.setZero();
      integral_grad_Jtprogress_to_ct.setZero();
      integral_grad_Jprogress_to_ct.setZero();
      integral_grad_Jtpositive_to_ct.setZero();
      integral_grad_Jtbound_to_ct.setZero();

      step_vaj = (total_t_ / N) / config_.K;
      s1 = 0.0;

      for (int j = 0; j <= config_.K; ++j, s1 += step_vaj)
      {
        //   rho = exp2(2 * (j + i * config_.K) / (N * config_.K));
        rho = 1;
        idx = i * config_.K + j;

        s2 = s1 * s1;
        s3 = s2 * s1;
        s4 = s2 * s2;
        s5 = s4 * s1;
        s6 = s4 * s2;
        s7 = s4 * s3;
        beta0 << 1.0, s1, s2, s3, s4, s5;
        beta1 << 0.0, 1.0, 2.0 * s1, 3.0 * s2, 4.0 * s3, 5.0 * s4;
        beta2 << 0.0, 0.0, 2.0, 6.0 * s1, 12.0 * s2, 20.0 * s3;
        beta3 << 0.0, 0.0, 0.0, 6.0, 24.0 * s1, 60.0 * s2;

        snap_beta0 << 1.0, s1, s2, s3, s4, s5, s6, s7;
        snap_beta1 << 0.0, 1.0, 2.0 * s1, 3.0 * s2, 4.0 * s3, 5.0 * s4, 6.0 * s5, 7.0 * s6;
        snap_beta2 << 0.0, 0.0, 2.0, 6.0 * s1, 12.0 * s2, 20.0 * s3, 30.0 * s4, 42.0 * s5;
        snap_beta3 << 0.0, 0.0, 0.0, 6.0, 24.0 * s1, 60.0 * s2, 120.0 * s3, 210.0 * s4;
        snap_beta4 << 0.0, 0.0, 0.0, 0.0, 24.0, 120.0 * s1, 360.0 * s2, 840.0 * s3;

        pos = c.transpose() * snap_beta0;
        vel = c.transpose() * snap_beta1;
        acc = c.transpose() * snap_beta2;
        jrk = c.transpose() * snap_beta3;

        t = ct.transpose() * beta0;
        t_dot = ct.transpose() * beta1;
        t_ddot = ct.transpose() * beta2;
        t_dddot = ct.transpose() * beta3;

        if (ref_traj_duration < t)
        {
          pos_ref = ref_traj_.getPos(ref_traj_duration);
          vel_ref = ref_traj_.getVel(ref_traj_duration);
        }
        else
        {
          pos_ref = ref_traj_.getPos(t);
          vel_ref = ref_traj_.getVel(t);
        }
        omega = (j == 0 || j == config_.K) ? 0.5 : 1.0;

        // cost construct
        t_progress_cost_ -= rho * step_vaj * omega * t_dot;
        integral_grad_Jtprogress_to_ct -= rho * omega * beta1;

        if (ProgressCostGrad(pos, pos_ref, grad_tmp, grad_ref_tmp, cost_tmp))
        {
          integral_grad_Jprogress_to_c += rho * omega * snap_beta0 * grad_tmp.transpose();
          integral_grad_Jprogress_to_ct += rho * omega * beta0 * vel_ref.transpose() * grad_ref_tmp;
          pos_progress_cost_ += rho * step_vaj * omega * cost_tmp;
        }

        // penalty construct
        // max limitation for pos traj
        if (NormalMaxConstriant(vel, grad_tmp, cost_tmp, config_.max_vel))
        {
          integral_grad_Jv_to_c += omega * snap_beta1 * grad_tmp.transpose();
          vel_pnt_ += step_vaj * omega * cost_tmp;
        }
        if (NormalMaxConstriant(acc, grad_tmp, cost_tmp, config_.max_acc))
        {
          integral_grad_Ja_to_c += omega * snap_beta2 * grad_tmp.transpose();
          acc_pnt_ += step_vaj * omega * cost_tmp;
        }
        if (NormalMaxConstriant(jrk, grad_tmp, cost_tmp, config_.max_jrk))
        {
          integral_grad_Jj_to_c += omega * snap_beta3 * grad_tmp.transpose();
          jrk_pnt_ += step_vaj * omega * cost_tmp;
        }

        //   max limitation for t traj
        if (NormalMaxConstriant1D(t_dot, tgrad_tmp, tcost_tmp, config_.T_max_vel))
        {
          integral_grad_Jtv_to_ct += omega * beta1 * tgrad_tmp;
          t_vel_pnt_ += step_vaj * omega * tcost_tmp;
        }
        if (NormalMaxConstriant1D(t_ddot, tgrad_tmp, tcost_tmp, config_.T_max_acc))
        {
          integral_grad_Jta_to_ct += omega * beta2 * tgrad_tmp;
          t_acc_pnt_ += step_vaj * omega * tcost_tmp;
        }

        if (NormalMaxConstriant1D(t_dddot, tgrad_tmp, tcost_tmp, config_.T_max_jrk))
        {
          integral_grad_Jtj_to_ct += omega * beta3 * tgrad_tmp;
          t_jrk_pnt_ += step_vaj * omega * tcost_tmp;
        }
        if (idx < (N - 1) * config_.K)
        {
          if (safeHpolyCostGrad(pos, grad_tmp, cost_tmp, CorridorPtr_->getHpoly(ids_[idx])))
          {
            integral_grad_JCone_to_c += omega * snap_beta0 * grad_tmp.transpose();
            Cone_pnt_ += step_vaj * omega * cost_tmp;
          }
        }

        if (PostiveTimeConstraint(t, tgrad_tmp, tcost_tmp))
        {
          integral_grad_Jtpositive_to_ct += omega * beta0 * tgrad_tmp;
          t_positive_pnt_ += step_vaj * omega * tcost_tmp;
        }

        // penalty to the time exceeding the time_horizon
        if (TimeboundCostGrad(t, tgrad_tmp, tcost_tmp))
        {
          integral_grad_Jtbound_to_ct += omega * beta0 * tgrad_tmp;
          t_bound_pnt_ += step_vaj * omega * tcost_tmp;
        }
      }

      gdC.block<8, 3>(i * 8, 0) +=
          step_vaj * (integral_grad_Jv_to_c * config_.pnlV +
                      integral_grad_Ja_to_c * config_.pnlA +
                      integral_grad_Jj_to_c * config_.pnlJ +
                      integral_grad_Jprogress_to_c * config_.rhoTrajProgress +
                      integral_grad_JCone_to_c * config_.pnlCone);

      gdCt.block<6, 1>(i * 6, 0) +=
          step_vaj * (integral_grad_Jtv_to_ct * config_.pnltV +
                      integral_grad_Jta_to_ct * config_.pnltA +
                      integral_grad_Jtj_to_ct * config_.pnltJ +
                      // integral_grad_Jtprogress_to_ct * config_.rhoTProgress +
                      integral_grad_Jtprogress_to_ct * rhoTprogress_ +
                      integral_grad_Jprogress_to_ct * config_.rhoTrajProgress +
                      integral_grad_Jtpositive_to_ct * config_.rhoTPositive +
                      integral_grad_Jtbound_to_ct * config_.rhoTup);
    }

    // add weights for penalty and cost
    //   t_progress_cost_ *= config_.rhoTProgress;
    t_progress_cost_ *= rhoTprogress_;

    pos_progress_cost_ *= config_.rhoTrajProgress;

    vel_pnt_ *= config_.pnlV;
    acc_pnt_ *= config_.pnlA;
    jrk_pnt_ *= config_.pnlJ;
    Cone_pnt_ *= config_.pnlCone;

    t_vel_pnt_ *= config_.pnltV;
    t_acc_pnt_ *= config_.pnltA;
    t_jrk_pnt_ *= config_.pnltJ;
    t_positive_pnt_ *= config_.rhoTPositive;
    t_bound_pnt_ *= config_.rhoTup;

    // std::cout << "vel_pnt : " << vel_pnt_ << "\n"
    //           << "acc_pnt : " << acc_pnt_ << "\n"
    //           << "jrk_pnt : " << jrk_pnt_ << "\n"
    //           << "t_vel_pnt : " << t_vel_pnt_ << "\n"
    //           << "t_acc_pnt : " << t_acc_pnt_ << "\n"
    //           << "t_jrk_pnt : " << t_jrk_pnt_ << "\n"
    //           << "t_progress_cost : " << t_progress_cost_ << "\n"
    //           << "pos_progress_cost : " << pos_progress_cost_ << "\n"
    //           << "t_positive_pnt_ : " << t_positive_pnt_ << "\n"
    //           << "Cone_pnt : " << Cone_pnt_ << std::endl;

    double accumulated_cost = vel_pnt_ + acc_pnt_ + jrk_pnt_ + Cone_pnt_ + t_progress_cost_ + pos_progress_cost_ + t_vel_pnt_ + t_acc_pnt_ + t_jrk_pnt_ + t_positive_pnt_ + t_bound_pnt_;

    return accumulated_cost;
  }

  bool UniformMincoProgress::TimeboundCostGrad(const double &t,
                                               double &grad,
                                               double &cost)
  {
    bool ret(false);
    cost = 0.0;
    grad = 0.0;
    double up_pnt = t - ini_T_[0] - total_t_;
    double lower_pnt = ini_T_[0] - t - total_t_;
    if (up_pnt > DBL_EPSILON)
    {
      double violaPena(up_pnt), violaPenaGrad(1.0);
      positiveSmoothedL1(up_pnt, config_.smoothEps, violaPena, violaPenaGrad);
      grad = violaPenaGrad;
      cost = violaPena;
      ret = true;
      return ret;
    }
    if (lower_pnt > DBL_EPSILON)
    {
      double violaPena(lower_pnt), violaPenaGrad(1.0);
      positiveSmoothedL1(lower_pnt, config_.smoothEps, violaPena, violaPenaGrad);
      grad = -violaPenaGrad;
      cost = violaPena;
      ret = true;
      return ret;
    }
    return ret;
  }

  bool UniformMincoProgress::safeHpolyCostGrad(const Eigen::Vector3d &pos,
                                               Eigen::Vector3d &grad,
                                               double &cost,
                                               Eigen::MatrixXd hpolys)
  {
    bool ret(false);
    cost = 0.0;
    grad.setZero();
    for (int i = 0; i < hpolys.cols(); i++)
    {
      Eigen::Vector3d p_, n_;
      p_ = hpolys.col(i).tail(3);
      n_ = hpolys.col(i).head(3);
      double pen = (pos - p_).dot(n_) + config_.corridor_safe_dist;
      if (pen > DBL_EPSILON)
      {
        double violaPena, violaPenaGrad;
        positiveSmoothedL1(pen, config_.smoothEps, violaPena, violaPenaGrad);
        grad += violaPenaGrad * n_;
        cost += violaPena;
        ret = true;
      }
    }

    return ret;
  }

  bool UniformMincoProgress::NormalMaxConstriant(const Eigen::Vector3d &derivative,
                                                 Eigen::Vector3d &grad,
                                                 double &cost,
                                                 const double max_limit)
  {
    bool ret(false);
    cost = 0.0;
    grad.setZero();
    double pen = derivative.squaredNorm() - max_limit * max_limit;
    if (pen > DBL_EPSILON)
    {
      double violaPena(pen), violaPenaGrad(1.0);
      positiveSmoothedL1(pen, config_.smoothEps, violaPena, violaPenaGrad);
      grad = 2.0 * derivative * violaPenaGrad;
      cost += violaPena;
      ret = true;
    }

    return ret;
  }

  bool UniformMincoProgress::NormalMaxConstriant1D(const double &derivative,
                                                   double &grad,
                                                   double &cost,
                                                   const double max_limit)
  {
    bool ret(false);
    cost = 0.0;
    grad = 0.0;
    double pen = derivative * derivative - max_limit * max_limit;
    if (pen > DBL_EPSILON)
    {
      double violaPena(pen), violaPenaGrad(1.0);
      positiveSmoothedL1(pen, config_.smoothEps, violaPena, violaPenaGrad);
      grad = 2.0 * derivative * violaPenaGrad;
      cost += violaPena;
      ret = true;
    }

    return ret;
  }

  bool UniformMincoProgress::PostiveTimeConstraint(const double &t,
                                                   double &grad,
                                                   double &cost)
  {
    bool ret(false);
    cost = 0.0;
    grad = 0.0;
    double pen = -t;
    if (pen > DBL_EPSILON)
    {
      double violaPena(pen), violaPenaGrad(1.0);
      positiveSmoothedL1(pen, config_.smoothEps, violaPena, violaPenaGrad);
      grad -= violaPenaGrad;
      cost += violaPena;
      ret = true;
    }

    return ret;
  }

  bool UniformMincoProgress::ProgressCostGrad(const Eigen::Vector3d &pos,
                                              const Eigen::Vector3d &ref_p,
                                              Eigen::Vector3d &gradp,
                                              Eigen::Vector3d &gradrefp,
                                              double &cost)
  {
    bool ret(false);
    cost = 0.0;
    gradp.setZero();
    gradrefp.setZero();
    // cost = (pos - ref_p).squaredNorm();
    double pen = (pos - ref_p).squaredNorm();
    if (pen > DBL_EPSILON)
    {
      cost += pen;
      gradp = 2 * (pos - ref_p);
      gradrefp = 2 * (ref_p - pos);
      ret = true;
    }
    return ret;
  }

  /*
  * For penalty of higher derivatives like vel\acc\jrk\...
  f = max(v^2 - v_max^2, 0)
  */
  bool UniformMincoProgress::highDerivativeCostGrad_vel(const Eigen::Vector3d &derivative,
                                                        Eigen::Vector3d &grad,
                                                        double &cost)
  {
    bool ret(false);
    cost = 0.0;
    grad.setZero();
    double sq_v_diff = derivative.squaredNorm() - squared_vel_limit_;
    if (sq_v_diff > DBL_EPSILON)
    {
      double violaPena(sq_v_diff), violaPenaGrad(1.0);
      positiveSmoothedL1(sq_v_diff, config_.smoothEps, violaPena, violaPenaGrad);
      grad = 2.0 * derivative * violaPenaGrad;
      cost += violaPena;
      ret = true;
    }

    return ret;
  }

  bool UniformMincoProgress::highDerivativeCostGrad_acc(const Eigen::Vector3d &derivative,
                                                        Eigen::Vector3d &grad,
                                                        double &cost)
  {
    bool ret(false);
    cost = 0.0;
    grad.setZero();
    // acceleration
    double sq_a_diff = derivative.squaredNorm() - squared_acc_limit_;
    if (sq_a_diff > DBL_EPSILON)
    {
      double violaPena(sq_a_diff), violaPenaGrad(1.0);
      positiveSmoothedL1(sq_a_diff, config_.smoothEps, violaPena, violaPenaGrad);
      grad = 2.0 * derivative * violaPenaGrad;
      cost += violaPena;
      ret = true;
    }

    return ret;
  }

  bool UniformMincoProgress::highDerivativeCostGrad_jrk(const Eigen::Vector3d &derivative,
                                                        Eigen::Vector3d &grad,
                                                        double &cost)
  {
    bool ret(false);
    cost = 0.0;
    grad.setZero();
    // jerk
    double sq_j_diff = derivative.squaredNorm() - squared_jrk_limit_;
    if (sq_j_diff > DBL_EPSILON)
    {
      double violaPena(sq_j_diff), violaPenaGrad(1.0);
      positiveSmoothedL1(sq_j_diff, config_.smoothEps, violaPena, violaPenaGrad);
      grad = 2.0 * derivative * violaPenaGrad;
      cost += violaPena;
      ret = true;
    }

    return ret;
  }

  // SECTION object function
  double UniformMincoProgress::objectiveFunc(void *ptrObj, const Eigen::VectorXd &x, Eigen::VectorXd &grad)
  {
    UniformMincoProgress *obj = reinterpret_cast<UniformMincoProgress *>(ptrObj);
    obj->iteration_++;
    Eigen::Map<const Eigen::MatrixXd> P1(x.data(), obj->waypoints_dim_, obj->N_1_ - 1);
    Eigen::Map<const Eigen::VectorXd> TP1(x.data() + obj->dim_p_1_, obj->dim_t_p_);
    Eigen::Map<const Eigen::MatrixXd> end_state(x.data() + obj->dim_p_1_ + obj->dim_t_p_, obj->waypoints_dim_, obj->state_order_);
    Eigen::Map<const Eigen::VectorXd> t_end_state(x.data() + obj->dim_p_1_ + obj->dim_t_p_ + obj->dim_end_state_, obj->dim_tend_state_);
    Eigen::Map<Eigen::MatrixXd> grad_p_1(grad.data(), obj->waypoints_dim_, obj->N_1_ - 1);
    Eigen::Map<Eigen::VectorXd> grad_t_p_1(grad.data() + obj->dim_p_1_, obj->dim_t_p_);
    Eigen::Map<Eigen::MatrixXd> grad_end_state(grad.data() + obj->dim_p_1_ + obj->dim_t_p_, obj->waypoints_dim_, obj->state_order_);
    Eigen::Map<Eigen::VectorXd> grad_t_end_state(grad.data() + obj->dim_p_1_ + obj->dim_t_p_ + obj->dim_end_state_, obj->dim_tend_state_);

    obj->snapOpt_1_.resetTailCon(end_state);
    obj->snapOpt_1_.generate(P1, obj->total_t_);

    obj->Timeopt_1_.resetTailCon(t_end_state);
    obj->Timeopt_1_.generate(TP1, obj->total_t_);

    double cost = 0.0;
    double energy_cost1, energy_cost2;
    // energy_cost1 = obj->jerkOpt_1_.getEnergy();

    energy_cost1 = obj->snapOpt_1_.getEnergy();
    energy_cost2 = obj->Timeopt_1_.getEnergy();
    // double energy_cost3 = obj->Timeopt_2_.getEnergy();
    cost += energy_cost1 + energy_cost2;

    obj->snapOpt_1_.getEnergyPartialGradByCoeffs(obj->gdC_1_);
    obj->Timeopt_1_.getEnergyPartialGradByCoeffs(obj->gdC_t_);
    double time_itl_cost_1 = obj->calTimeIntPenalty(obj->N_1_, obj->gdC_1_, obj->gdC_t_);
    cost += time_itl_cost_1;

    double assist_cost = obj->calFixedTimeIntPenalty(obj->N_1_, obj->gdC_1_);
    cost += assist_cost;

    obj->snapOpt_1_.propogateGrad(obj->gdC_1_, obj->gdP_1_);
    obj->Timeopt_1_.propogateGrad(obj->gdC_t_, obj->gdP_t_);

    Eigen::MatrixXd tmp_mid_state_grad;
    obj->snapOpt_1_.getPartailGradsByEndState(tmp_mid_state_grad);
    grad_end_state = tmp_mid_state_grad.transpose();

    Eigen::VectorXd tmp_t_state_grad;
    obj->Timeopt_1_.getPartailGradsByEndState(tmp_t_state_grad);
    grad_t_end_state = tmp_t_state_grad;

    grad_p_1 = obj->gdP_1_;
    grad_t_p_1 = obj->gdP_t_;

    return cost;
  }

  double UniformMincoProgress::calFixedTimeIntPenalty(const int &N,
                                                      Eigen::MatrixX3d &gdC)
  {
    Eigen::Vector3d pos_grad;
    Eigen::Vector4d grad_q_atti;
    Eigen::Vector3d grad_body_rate;
    Eigen::Vector3d grad_Jfm_to_vel;
    Eigen::Vector3d grad_Jfm_to_acc;
    Eigen::Vector3d grad_Jfm_to_jrk;

    Eigen::Vector3d pos, vel, acc, jrk;
    Eigen::Vector3d grad_tmp;
    double cost_tmp;
    static Eigen::Matrix<double, 8, 1> beta0, beta1, beta2, beta3, beta4;
    double s1, s2, s3, s4, s5, s6, s7;
    double omega;
    Eigen::Matrix<double, 8, 3> grad_Jp_to_c, grad_Jv_to_c, grad_Ja_to_c, grad_Jj_to_c;
    Eigen::Matrix<double, 8, 3> integral_grad_Jp_to_c, integral_grad_Jv_to_c, integral_grad_Ja_to_c, integral_grad_Jj_to_c, integral_grad_Jfm_to_c;
    double grad_Jp_to_T, grad_Jv_to_T, grad_Ja_to_T;
    double integral_grad_Jp_to_T, integral_grad_Jv_to_T, integral_grad_Ja_to_T;
    double integral_cost_p, integral_cost_v, integral_cost_a;

    const Eigen::MatrixX3d &coeffs = snapOpt_1_.getCoeffs();

    int K_fix = K_fix_;
    double step_fix = step_fix_;

    int IdxPiece = 0;
    double t = 0.0;
    double t_pre = 0.0;
    double dT = total_t_ / N;
    integral_cost_a = 0;
    integral_cost_v = 0;
    integral_cost_p = 0;

    for (int j = 0; j <= K_fix; ++j)
    {
      integral_grad_Ja_to_T = 0;
      integral_grad_Jv_to_T = 0;
      integral_grad_Jp_to_T = 0;

      integral_grad_Jp_to_c.setZero();
      integral_grad_Jv_to_c.setZero();
      integral_grad_Ja_to_c.setZero();

      double rho = 1.0;
      while (t - t_pre > dT)
      {
        t_pre += dT;
        IdxPiece++;
      }
      const auto &c = coeffs.block<8, 3>(IdxPiece * 8, 0);
      s1 = t - t_pre;
      s2 = s1 * s1;
      s3 = s2 * s1;
      s4 = s2 * s2;
      s5 = s4 * s1;
      s6 = s4 * s2;
      s7 = s4 * s3;

      beta0 << 1.0, s1, s2, s3, s4, s5, s6, s7;
      beta1 << 0.0, 1.0, 2.0 * s1, 3.0 * s2, 4.0 * s3, 5.0 * s4, 6.0 * s5, 7.0 * s6;
      beta2 << 0.0, 0.0, 2.0, 6.0 * s1, 12.0 * s2, 20.0 * s3, 30.0 * s4, 42.0 * s5;
      beta3 << 0.0, 0.0, 0.0, 6.0, 24.0 * s1, 60.0 * s2, 120.0 * s3, 210.0 * s4;
      beta4 << 0.0, 0.0, 0.0, 0.0, 24.0, 120.0 * s1, 360.0 * s2, 840.0 * s3;
      pos = c.transpose() * beta0;
      vel = c.transpose() * beta1;
      acc = c.transpose() * beta2;
      jrk = c.transpose() * beta3;

      Eigen::Vector3d acc_des = Prim.getAcc(t);
      Eigen::Vector3d vel_des = Prim.getVel(t);
      Eigen::Vector3d pos_des = Prim.getPos(t);
      omega = (j == 0 || j == K_fix) ? 0.5 : 1.0;

      if (d_set_[j] == -1)
      {
        t += step_fix;
        break;
      }
      else
      {
        if (d_set_[j] < config_.corridor_edge_release_dist)
          rho *= exp2(-8.0 * (config_.corridor_edge_release_dist - d_set_[j]));
      }

      if (UserIntentionPenalty(acc, acc_des, grad_tmp, cost_tmp))
      {
        grad_tmp = config_.rhoAccCmd * grad_tmp;
        cost_tmp = config_.rhoAccCmd * cost_tmp;
        grad_Ja_to_c = beta2 * grad_tmp.transpose();
        grad_Ja_to_T = grad_tmp.dot(jrk);
        integral_grad_Ja_to_c += omega * step_fix * grad_Ja_to_c;
        if (IdxPiece > 0)
        {
          integral_grad_Ja_to_T += -IdxPiece * omega * step_fix * grad_Ja_to_T;
        }
        integral_cost_a += rho * step_fix * omega * cost_tmp;
      }
      if (UserIntentionPenalty(vel, vel_des, grad_tmp, cost_tmp))
      {
        grad_tmp = config_.rhoVelCmd * grad_tmp;
        cost_tmp = config_.rhoVelCmd * cost_tmp;
        grad_Jv_to_c = beta1 * grad_tmp.transpose();
        grad_Jv_to_T = grad_tmp.dot(acc);
        integral_grad_Jv_to_c += omega * step_fix * grad_Jv_to_c;
        if (IdxPiece > 0)
        {
          integral_grad_Jv_to_T += -IdxPiece * omega * step_fix * grad_Jv_to_T;
        }
        integral_cost_v += rho * step_fix * omega * cost_tmp;
      }

      if (UserIntentionPenalty(pos, pos_des, grad_tmp, cost_tmp))
      {
        // std::cout<<"tragger!acc: "<<acc.transpose()<<std::endl;
        grad_tmp = config_.rhoPosCmd * grad_tmp;
        cost_tmp = config_.rhoPosCmd * cost_tmp;
        grad_Jp_to_c = beta0 * grad_tmp.transpose();
        grad_Jp_to_T = grad_tmp.dot(vel);
        integral_grad_Jp_to_c += omega * step_fix * grad_Jp_to_c;
        if (IdxPiece > 0)
        {
          integral_grad_Jp_to_T += -IdxPiece * omega * step_fix * grad_Jp_to_T;
        }
        integral_cost_p += rho * step_fix * omega * cost_tmp;
      }
      gdC.block<8, 3>(IdxPiece * 8, 0) += rho * (integral_grad_Ja_to_c + integral_grad_Jv_to_c + integral_grad_Jp_to_c);

      t += step_fix;
    }
    Rc_pos_pnt_ = integral_cost_p;
    Rc_vel_pnt_ = integral_cost_v;
    Rc_acc_pnt_ = integral_cost_a;

    return (integral_cost_a + integral_cost_v + integral_cost_p);
  }

  bool UniformMincoProgress::UserIntentionPenalty(const Eigen::Vector3d &acc,
                                                  Eigen::Vector3d &des_acc,
                                                  Eigen::Vector3d &grad,
                                                  double &cost)
  {
    bool ret(false);
    cost = 0.0;
    grad.setZero();
    for (int i = 0; i < 3; i++)
    {
      double pen = (acc[i] - des_acc[i]) * (acc[i] - des_acc[i]);
      if (pen > 0)
      {
        double violaPena, violaPenaGrad;
        positiveSmoothedL1(pen, config_.smoothEps, violaPena, violaPenaGrad);
        grad[i] = 2 * (acc[i] - des_acc[i]) * violaPenaGrad;
        cost += violaPena;
        ret = true;
      }
    }

    return ret;
  }

  inline int UniformMincoProgress::Process(void *ptrObj, const double *x,
                                           const double *grad, const double fx,
                                           const double xnorm, const double gnorm,
                                           const double step, int n, int k, int ls)
  {
    UniformMincoProgress *obj = reinterpret_cast<UniformMincoProgress *>(ptrObj);

    if (obj->config_.debug == false)
    {
      return 0;
    }

    ros::Duration(0.01).sleep();

    return 0;
  }

  inline void UniformMincoProgress::forwardT(const Eigen::Ref<const Eigen::VectorXd> &t,
                                             Eigen::Ref<Eigen::VectorXd> vecT)
  {
    int M = t.size();
    for (int i = 0; i < M; ++i)
    {
      vecT(i) = expC2(t(i));
    }
    return;
  }

  inline void UniformMincoProgress::backwardT(
      const Eigen::Ref<const Eigen::VectorXd> &vecT,
      Eigen::Ref<Eigen::VectorXd> t)
  {
    int M = t.size();
    for (int i = 0; i < M; ++i)
    {
      t(i) = logC2(vecT(i));
    }
    return;
  }

  inline void UniformMincoProgress::addLayerTGrad(
      const Eigen::Ref<const Eigen::VectorXd> &t,
      const Eigen::Ref<const Eigen::VectorXd> &gradT,
      Eigen::Ref<Eigen::VectorXd> gradt)
  {
    int dim_t = t.size();
    for (int i = 0; i < dim_t; ++i)
    {
      gradt[i] = gradT[i] * d_quasi_exp(t[i]);
    }
  }

  inline void UniformMincoProgress::addLayerPGrad(
      const Eigen::Ref<const Eigen::MatrixXd> &gradInPs,
      Eigen::Ref<Eigen::MatrixXd> grad)
  {
    grad = gradInPs;
    return;
  }

  // this is a smoothed C2 exact penalty
  // inputs x >0, output exact penalty and penalty derivates
  inline void UniformMincoProgress::positiveSmoothedL1(const double &x, const double &seps,
                                                       double &f, double &df)
  {
    static const double pe = seps;
    static const double half = 0.5 * pe;
    static const double f3c = 1.0 / (pe * pe);
    static const double f4c = -0.5 * f3c / pe;
    static const double d2c = 3.0 * f3c;
    static const double d3c = 4.0 * f4c;

    if (x < pe)
    {
      f = (f4c * x + f3c) * x * x * x;
      df = (d3c * x + d2c) * x * x;
    }
    else
    {
      f = x - half;
      df = 1.0;
    }

    return;
  }

  inline double UniformMincoProgress::d_quasi_exp(const double &x)
  {
    if (x > 0.0)
    {
      return x + 1.0;
    }
    else
    {
      double denSqrt = (0.5 * x - 1.0) * x + 1.0;
      return (1.0 - x) / (denSqrt * denSqrt);
    }
  }

  double UniformMincoProgress::expC2(double t)
  {
    return t > 0.0 ? ((0.5 * t + 1.0) * t + 1.0)
                   : 1.0 / ((0.5 * t - 1.0) * t + 1.0);
  }

  double UniformMincoProgress::logC2(double T)
  {
    return T > 1.0 ? (sqrt(2.0 * T - 1.0) - 1.0) : (1.0 - sqrt(2.0 / T - 1.0));
  }

} // namespace traj_opt
