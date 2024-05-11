#include <traj_opt/uniform_minco_emergency.h>

#include <traj_opt/chalk.hpp>
#include <optimization_opt/lbfgs_raw.hpp>

namespace traj_opt
{
  UniformMincoStop::UniformMincoStop(ros::NodeHandle &nh, Config &conf)
      : nh_(nh), config_(conf), gdT_(2), propagated_gdT_(2)
  {
    squared_vel_limit_ = config_.max_vel * config_.max_vel;
    squared_acc_limit_ = config_.max_acc * config_.max_acc;
    squared_jrk_limit_ = config_.max_jrk * config_.max_jrk;
  }

  bool UniformMincoStop::generate_traj(const Eigen::MatrixXd &iniState,
                                       const Eigen::MatrixXd &finState,
                                       const Eigen::MatrixXd &waypoints_1,
                                       const double &dt_1,
                                       Trajectory &traj)
  {
    ini_state_ = iniState;
    end_state_ = finState;

    waypoints_1_ = waypoints_1;

    waypoints_dim_ = waypoints_1.rows();
    dim_p_1_ = waypoints_dim_ * waypoints_1.cols();

    N_1_ = waypoints_1.cols() + 1;

    state_order_ = finState.cols();
    dim_end_state_ = waypoints_dim_ * state_order_;

    dim_t_ = 1;

    Eigen::VectorXd total_T;
    total_T.resize(dim_t_);
    total_T << dt_1 * N_1_;

    roundingState(ini_state_);
    roundingState(end_state_);
    setBoundConds(ini_state_, end_state_);

    x_.resize(dim_p_1_ + dim_t_ + dim_end_state_);
    Eigen::Map<Eigen::MatrixXd> P1(x_.data(), waypoints_dim_, N_1_ - 1);
    Eigen::Map<Eigen::VectorXd> total_t(x_.data() + dim_p_1_, dim_t_);
    Eigen::Map<Eigen::MatrixXd> end_state(x_.data() + dim_p_1_ + dim_t_, waypoints_dim_, state_order_);

    P1 = waypoints_1;
    backwardT(total_T, total_t);
    end_state = end_state_;

    cost_function_cb_times_ = 0;

    jerkOpt_1_.generate(P1, total_T[0]);

    int opt_ret = optimize();

    forwardT(total_t, total_T);

    if (opt_ret < 0)
    {
      return false;
    }

    jerkOpt_1_.generate(P1, total_T[0]);

    traj = jerkOpt_1_.getTraj();

    return true;
  }

  int UniformMincoStop::optimize(const double &delta)
  {
    // Setup for L-BFGS solver
    lbfgs::lbfgs_parameter_t lbfgs_params;
    lbfgs_params.mem_size = 32;
    lbfgs_params.past = 3;
    lbfgs_params.g_epsilon = 0.0;
    lbfgs_params.min_step = 1e-32;
    lbfgs_params.delta = 1.0e-6;
    double minObjective;
    auto ret = lbfgs::lbfgs_optimize(
        x_, minObjective, UniformMincoStop::objectiveFunc, nullptr,
        nullptr, this, lbfgs_params);
    std::cout << "\033[32m"
              << "ret: " << ret << ", minObjective: " << minObjective << "\033[0m"
              << std::endl;
    return ret;
  }

  inline void UniformMincoStop::roundingState(Eigen::MatrixXd &state)
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
  }

  inline void UniformMincoStop::setBoundConds(const Eigen::MatrixXd &iniState,
                                              const Eigen::MatrixXd &finState)
  {
    jerkOpt_1_.reset(iniState, finState, N_1_);
  }

  double UniformMincoStop::calGatePenalty(const Eigen::MatrixXd &state, Eigen::Ref<Eigen::MatrixXd> grad)
  {
    double cost = 0.0;
    double violation = (state.col(0) - gate_pos_).squaredNorm();
    cost = config_.pnlGate * violation;
    grad.col(0) += config_.pnlGate * 2.0 * (state.col(0) - gate_pos_);
    return cost;
  }

  double UniformMincoStop::calTimeIntPenalty(const int &N,
                                             Eigen::MatrixX3d &gdC,
                                             double &gdT,
                                             double total_T,
                                             const Eigen::MatrixX3d &coeffs)
  {
    Eigen::Vector3d pos, vel, acc, jrk, snp;
    Eigen::Vector3d grad_tmp;
    double cost_tmp;
    static Eigen::Matrix<double, 6, 1> beta0, beta1, beta2, beta3, beta4;
    double s1, s2, s3, s4, s5;
    double step_vaj, omega, alpha;
    Eigen::Matrix<double, 6, 3> grad_Jv_to_c, grad_Ja_to_c, grad_Jj_to_c, grad_JCone_to_c;
    Eigen::Matrix<double, 6, 3> integral_grad_Jv_to_c, integral_grad_Ja_to_c, integral_grad_JCone_to_c, integral_grad_Jj_to_c;
    double grad_Jv_to_T, grad_Ja_to_T, grad_Jj_to_T, grad_JCone_to_T;
    double integral_grad_Jv_to_T, integral_grad_Ja_to_T, integral_grad_Jj_to_T, integral_grad_JCone_to_T;
    double integral_cost_v, integral_cost_a, integral_cost_j, integral_cost_Cone;

    vel_pnt_ = 0.0;
    acc_pnt_ = 0.0;
    jrk_pnt_ = 0.0;
    Cone_pnt_ = 0.0;
    double accumulated_dur(0.0);
    for (int i = 0; i < N; ++i)
    {
      const auto &c = coeffs.block<6, 3>(i * 6, 0);
      integral_grad_Jv_to_c.setZero();
      integral_grad_Ja_to_c.setZero();
      integral_grad_Jj_to_c.setZero();
      integral_grad_JCone_to_c.setZero();
      integral_grad_Jv_to_T = 0.0;
      integral_grad_Ja_to_T = 0.0;
      integral_grad_Jj_to_T = 0.0;
      integral_grad_JCone_to_T = 0.0;
      integral_cost_v = 0.0;
      integral_cost_a = 0.0;
      integral_cost_j = 0.0;
      integral_cost_Cone = 0.0;

      step_vaj = (total_T / N) / config_.K;
      s1 = 0.0;
      for (int j = 0; j <= config_.K; ++j, s1 += step_vaj)
      {
        s2 = s1 * s1;
        s3 = s2 * s1;
        s4 = s2 * s2;
        s5 = s4 * s1;
        beta0[0] = 1.0;
        beta0[1] = s1;
        beta0[2] = s2;
        beta0[3] = s3;
        beta0[4] = s4;
        beta0[5] = s5;
        beta1[0] = 0.0;
        beta1[1] = 1.0;
        beta1[2] = 2.0 * s1;
        beta1[3] = 3.0 * s2;
        beta1[4] = 4.0 * s3;
        beta1[5] = 5.0 * s4;
        beta2[0] = 0.0;
        beta2[1] = 0.0;
        beta2[2] = 2.0;
        beta2[3] = 6.0 * s1;
        beta2[4] = 12.0 * s2;
        beta2[5] = 20.0 * s3;
        beta3[0] = 0.0;
        beta3[1] = 0.0;
        beta3[2] = 0.0;
        beta3[3] = 6.0;
        beta3[4] = 24.0 * s1;
        beta3[5] = 60.0 * s2;
        beta4[0] = 0.0;
        beta4[1] = 0.0;
        beta4[2] = 0.0;
        beta4[3] = 0.0;
        beta4[4] = 24.0;
        beta4[5] = 120.0 * s1;
        pos = c.transpose() * beta0;
        vel = c.transpose() * beta1;
        acc = c.transpose() * beta2;
        jrk = c.transpose() * beta3;
        snp = c.transpose() * beta4;
        alpha = (double)j / config_.K;
        omega = (j == 0 || j == config_.K) ? 0.5 : 1.0;
        if (highDerivativeCostGrad_vel(vel, grad_tmp, cost_tmp))
        {
          grad_Jv_to_c = beta1 * grad_tmp.transpose();
          grad_Jv_to_T = alpha * grad_tmp.dot(acc);
          integral_grad_Jv_to_c += omega * grad_Jv_to_c;
          integral_grad_Jv_to_T +=
              omega * (cost_tmp / config_.K + step_vaj * grad_Jv_to_T);
          integral_cost_v += omega * cost_tmp;
        }
        if (highDerivativeCostGrad_acc(acc, grad_tmp, cost_tmp))
        {
          grad_Ja_to_c = beta2 * grad_tmp.transpose();
          grad_Ja_to_T = alpha * grad_tmp.dot(jrk);
          integral_grad_Ja_to_c += omega * grad_Ja_to_c;
          integral_grad_Ja_to_T +=
              omega * (cost_tmp / config_.K + step_vaj * grad_Ja_to_T);
          integral_cost_a += omega * cost_tmp;
        }
        if (highDerivativeCostGrad_jrk(jrk, grad_tmp, cost_tmp))
        {
          grad_Jj_to_c = beta3 * grad_tmp.transpose();
          grad_Jj_to_T = alpha * grad_tmp.dot(snp);
          integral_grad_Jj_to_c += omega * grad_Jj_to_c;
          integral_grad_Jj_to_T +=
              omega * (cost_tmp / config_.K + step_vaj * grad_Jj_to_T);
          integral_cost_j += omega * cost_tmp;
        }

        if (safeHpolyCostGrad(pos, grad_tmp, cost_tmp, Safe_box_))
        {
          grad_JCone_to_c = beta0 * grad_tmp.transpose();
          grad_JCone_to_T = alpha * grad_tmp.dot(vel);
          integral_grad_JCone_to_c += omega * grad_JCone_to_c;
          integral_grad_JCone_to_T +=
              omega * (cost_tmp / config_.K + step_vaj * grad_JCone_to_T);
          integral_cost_Cone += omega * cost_tmp;
        }
      }
      vel_pnt_ += step_vaj * integral_cost_v;
      acc_pnt_ += step_vaj * integral_cost_a;
      jrk_pnt_ += step_vaj * integral_cost_j;
      Cone_pnt_ += step_vaj * integral_cost_Cone;

      gdC.block<6, 3>(i * 6, 0) +=
          step_vaj * (integral_grad_Jv_to_c * config_.pnlV +
                      integral_grad_Ja_to_c * config_.pnlA +
                      integral_grad_Jj_to_c * config_.pnlJ +
                      integral_grad_JCone_to_c * config_.pnlCone);
      gdT += (integral_grad_Jv_to_T * config_.pnlV +
              integral_grad_Ja_to_T * config_.pnlA +
              integral_grad_Jj_to_T * config_.pnlJ +
              integral_grad_JCone_to_T * config_.pnlCone) /
             N;

      accumulated_dur += total_T / N;
    }
    vel_pnt_ *= config_.pnlV;
    acc_pnt_ *= config_.pnlA;
    jrk_pnt_ *= config_.pnlJ;
    Cone_pnt_ *= config_.pnlCone;

    // std::cout << "vel_pnt : " << vel_pnt_ << "\n"
    //           << "acc_pnt : " << acc_pnt_ << "\n"
    //           << "jrk_pnt : " << jrk_pnt_ << std::endl;

    return (vel_pnt_ + acc_pnt_ + jrk_pnt_ + Cone_pnt_);
  }

  bool UniformMincoStop::safeHpolyCostGrad(const Eigen::Vector3d &pos,
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

  /*
  * For penalty of higher derivatives like vel\acc\jrk\...
  f = max(v^2 - v_max^2, 0)
  */
  bool UniformMincoStop::highDerivativeCostGrad_vel(const Eigen::Vector3d &derivative,
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

  bool UniformMincoStop::highDerivativeCostGrad_acc(const Eigen::Vector3d &derivative,
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

  bool UniformMincoStop::highDerivativeCostGrad_jrk(const Eigen::Vector3d &derivative,
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
  double UniformMincoStop::objectiveFunc(void *ptrObj, const Eigen::VectorXd &x, Eigen::VectorXd &grad)
  {
    UniformMincoStop *obj = reinterpret_cast<UniformMincoStop *>(ptrObj);

    Eigen::Map<const Eigen::MatrixXd> P1(x.data(), obj->waypoints_dim_, obj->N_1_ - 1);
    Eigen::Map<const Eigen::VectorXd> total_t(x.data() + obj->dim_p_1_, obj->dim_t_);
    Eigen::Map<const Eigen::MatrixXd> end_state(x.data() + obj->dim_p_1_ + obj->dim_t_, obj->waypoints_dim_, obj->state_order_);
    Eigen::Map<Eigen::MatrixXd> grad_p_1(grad.data(), obj->waypoints_dim_, obj->N_1_ - 1);
    Eigen::Map<Eigen::VectorXd> grad_total_t(grad.data() + obj->dim_p_1_, obj->dim_t_);
    Eigen::Map<Eigen::MatrixXd> grad_end_state(grad.data() + obj->dim_p_1_ + obj->dim_t_, obj->waypoints_dim_, obj->state_order_);

    // std::cout<<"mid_state: "<< mid_state <<std::endl;

    Eigen::VectorXd total_T;
    total_T.resize(obj->dim_t_);
    forwardT(total_t, total_T);

    obj->jerkOpt_1_.resetTailCon(end_state);
    obj->jerkOpt_1_.generate(P1, total_T[0]);

    double cost = 0.0;
    cost += obj->jerkOpt_1_.getEnergy();

    // std::cout << "--------------------------------------------------" << std::endl;
    // std::cout << "energy_cost : " << cost << std::endl;

    obj->jerkOpt_1_.getEnergyPartialGradByCoeffs(obj->gdC_1_);
    obj->jerkOpt_1_.getEnergyPartialGradByTotalTime(obj->gdT_[0]);

    double time_itl_cost_1 = obj->calTimeIntPenalty(obj->N_1_, obj->gdC_1_, obj->gdT_[0], total_T[0], obj->jerkOpt_1_.getCoeffs());

    cost += time_itl_cost_1;

    obj->jerkOpt_1_.propogateGrad(obj->gdC_1_, obj->gdT_[0], obj->gdP_1_, obj->propagated_gdT_[0]);

    Eigen::MatrixXd tmp_mid_state_grad;
    tmp_mid_state_grad.resize(obj->waypoints_dim_, obj->state_order_);

    grad_end_state.setZero();

    obj->jerkOpt_1_.getPartailGradsByEndState(tmp_mid_state_grad);
    grad_end_state += tmp_mid_state_grad.transpose();

    // grad_end_state.row(2).setZero();
    grad_end_state.col(1).setZero();
    grad_end_state.col(2).setZero();

    // cost += obj->calGatePenalty(mid_state, grad_mid_state);

    cost += obj->config_.rhoT * total_T.sum();
    obj->propagated_gdT_.array() += obj->config_.rhoT;

    addLayerTGrad(total_t, obj->propagated_gdT_, grad_total_t);
    addLayerPGrad(obj->gdP_1_, grad_p_1);

    return cost;
  }

  inline int UniformMincoStop::Process(void *ptrObj, const double *x,
                                       const double *grad, const double fx,
                                       const double xnorm, const double gnorm,
                                       const double step, int n, int k, int ls)
  {
    UniformMincoStop *obj = reinterpret_cast<UniformMincoStop *>(ptrObj);

    if (obj->config_.debug == false)
    {
      return 0;
    }

    Eigen::Map<const Eigen::MatrixXd> P1(x, obj->waypoints_dim_, obj->N_1_ - 1);
    Eigen::Map<const Eigen::VectorXd> total_t(x + obj->dim_p_1_, obj->dim_t_);
    Eigen::Map<const Eigen::MatrixXd> end_state(x + obj->dim_p_1_ + obj->dim_t_, obj->waypoints_dim_, obj->state_order_);

    Eigen::VectorXd total_T;
    total_T.resize(obj->dim_t_);
    forwardT(total_t, total_T);

    obj->jerkOpt_1_.resetTailCon(end_state);

    obj->jerkOpt_1_.generate(P1, total_T[0]);

    Trajectory traj;
    traj = obj->jerkOpt_1_.getTraj();

    obj->visPtr_->visualize_traj(traj, "debug_traj");

    ros::Duration(0.01).sleep();

    return 0;
  }

  inline void UniformMincoStop::forwardT(const Eigen::Ref<const Eigen::VectorXd> &t,
                                         Eigen::Ref<Eigen::VectorXd> vecT)
  {
    int M = t.size();
    for (int i = 0; i < M; ++i)
    {
      vecT(i) = expC2(t(i));
    }
    return;
  }

  inline void UniformMincoStop::backwardT(
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

  inline void UniformMincoStop::addLayerTGrad(
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

  inline void UniformMincoStop::addLayerPGrad(
      const Eigen::Ref<const Eigen::MatrixXd> &gradInPs,
      Eigen::Ref<Eigen::MatrixXd> grad)
  {
    grad = gradInPs;
    return;
  }

  // this is a smoothed C2 exact penalty
  // inputs x >0, output exact penalty and penalty derivates
  inline void UniformMincoStop::positiveSmoothedL1(const double &x, const double &seps,
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
  inline double UniformMincoStop::d_quasi_exp(const double &x)
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
    // return 2.0 * x;
  }
  double UniformMincoStop::expC2(double t)
  {
    return t > 0.0 ? ((0.5 * t + 1.0) * t + 1.0)
                   : 1.0 / ((0.5 * t - 1.0) * t + 1.0);
    // return 0.5 + t * t;
  }
  double UniformMincoStop::logC2(double T)
  {
    return T > 1.0 ? (sqrt(2.0 * T - 1.0) - 1.0) : (1.0 - sqrt(2.0 / T - 1.0));
    // return sqrt(std::max(T - 0.5, 0.0));
  }

} // namespace traj_opt
