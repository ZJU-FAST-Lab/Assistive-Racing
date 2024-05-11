#include <traj_opt/global_traj_opt.h>

#include <traj_opt/chalk.hpp>
#include <optimization_opt/lbfgs_raw.hpp>

namespace traj_opt
{
  GlobTrajOpt::GlobTrajOpt(ros::NodeHandle &nh, Config &conf)
      : nh_(nh), config_(conf)
  {
    squared_vel_limit_ = config_.max_vel * config_.max_vel;
    squared_acc_limit_ = config_.max_acc * config_.max_acc;
    squared_jrk_limit_ = config_.max_jrk * config_.max_jrk;
    thrAccMinSqr = config_.min_thracc * config_.min_thracc;
    thrAccMaxSqr = config_.max_thracc * config_.max_thracc;
  }

  bool GlobTrajOpt::generate_traj(const Eigen::MatrixXd &iniState,
                                  const Eigen::MatrixXd &finState,
                                  const Eigen::MatrixXd &waypoints,
                                  const Eigen::VectorXd &Ts, Trajectory &traj)
  {
    ini_state_ = iniState;
    end_state_ = finState;
    iter_num = 0;
    double Ts_sum(Ts.sum());
    printf(CHALK_YELLOW("original Ts: %lf\n"), Ts_sum);

    // initialize waypoints
    waypoints_ = waypoints;

    // initialize time durations
    N_ = Ts.size();
    dim_t_ = N_;
    t_.resize(dim_t_);
    backwardT(Ts, t_);

    roundingState(ini_state_);
    roundingState(end_state_);
    setBoundConds(ini_state_, end_state_);
    // we optimize the end velocity and accleration
    x_.resize(dim_t_ + 6);
    Eigen::Map<Eigen::VectorXd> t(x_.data(), dim_t_);
    t = t_;
    memcpy(x_.data() + dim_t_, end_state_.col(1).data(), 3 * sizeof(x_[0]));     // end vel
    memcpy(x_.data() + dim_t_ + 3, end_state_.col(2).data(), 3 * sizeof(x_[0])); // end acc

    int opt_ret = optimize();
    if (opt_ret < 0)
    {
      return false;
    }

    Eigen::VectorXd T(N_);
    forwardT(t, T);

    printf(CHALK_GREEN("optimized Ts: %lf\n"), T.sum());

    jerkOpt_.generate(waypoints, T);
    traj = jerkOpt_.getTraj();
    return true;
  }

  double GlobTrajOpt::calFixedTimeIntPenalty()
  {

    Eigen::Vector3d pos_grad;
    Eigen::Vector3d grad_body_rate;
    Eigen::Vector3d grad_Jfm_to_vel;
    Eigen::Vector3d grad_Jfm_to_acc;
    Eigen::Vector3d grad_Jfm_to_jrk;

    Eigen::Vector3d pos, vel, acc, jrk, snp;
    Eigen::Vector3d grad_tmp;
    double cost_tmp;
    static Eigen::Matrix<double, 6, 1> beta0, beta1, beta2, beta3, beta4;
    double s1, s2, s3, s4, s5;
    double omega;
    Eigen::Matrix<double, 6, 3> grad_Jp_to_c, grad_Jv_to_c, grad_Ja_to_c, grad_Jj_to_c;
    Eigen::Matrix<double, 6, 3> integral_grad_Jp_to_c, integral_grad_Jv_to_c, integral_grad_Ja_to_c, integral_grad_Jj_to_c, integral_grad_Jfm_to_c;
    double grad_Ja_to_T;
    double integral_cost_a;

    int K_fix = 20;
    double step_fix = 0.01;

    int IdxPiece = 0;
    double t = 0.0;
    double t_pre = 0.0;
    integral_cost_a = 0.0;
    integral_grad_Ja_to_c.setZero();

    for (int j = 0; j <= K_fix; ++j)
    {

      while (t - t_pre > jerkOpt_.T1(IdxPiece))
      {
        t_pre += jerkOpt_.T1(IdxPiece);
        IdxPiece++;
      }

      const auto &c = jerkOpt_.b.block<6, 3>(IdxPiece * 6, 0);
      s1 = t - t_pre;
      s2 = s1 * s1;
      s3 = s2 * s1;
      s4 = s2 * s2;
      s5 = s4 * s1;
      beta0 << 1.0, s1, s2, s3, s4, s5;
      beta1 << 0.0, 1.0, 2.0 * s1, 3.0 * s2, 4.0 * s3, 5.0 * s4;
      beta2 << 0.0, 0.0, 2.0, 6.0 * s1, 12.0 * s2, 20.0 * s3;
      beta3 << 0.0, 0.0, 0.0, 6.0, 24.0 * s1, 60.0 * s2;
      beta4 << 0.0, 0.0, 0.0, 0.0, 24.0, 120.0 * s1;
      pos = c.transpose() * beta0;
      vel = c.transpose() * beta1;
      acc = c.transpose() * beta2;
      jrk = c.transpose() * beta3;
      snp = c.transpose() * beta4;

      Eigen::Vector3d acc_des(0, 0, 0);
      omega = (j == 0 || j == K_fix) ? 0.5 : 1.0;

      if (UserIntentionPenalty(acc, acc_des, grad_tmp, cost_tmp))
      {
        grad_Ja_to_c = beta2 * grad_tmp.transpose();
        grad_Ja_to_T = grad_tmp.dot(jrk);
        integral_grad_Ja_to_c += omega * step_fix * grad_Ja_to_c;
        if (IdxPiece > 0)
        {
          jerkOpt_.gdT.head(IdxPiece).array() += -omega * step_fix * grad_Ja_to_T;
        }
        integral_cost_a += step_fix * omega * cost_tmp;
      }
      jerkOpt_.gdC.block<6, 3>(IdxPiece * 6, 0) += integral_grad_Ja_to_c;

      t += step_fix;
    }
    return integral_cost_a;
  }

  bool GlobTrajOpt::UserIntentionPenalty(const Eigen::Vector3d &acc,
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
        grad[i] = rhoCmdAcc_ * 2 * (acc[i] - des_acc[i]) * violaPenaGrad;
        cost += rhoCmdAcc_ * violaPena;
        ret = true;
      }
    }

    return ret;
  }

  int GlobTrajOpt::optimize(const double &delta)
  {
    // Setup for L-BFGS solver
    lbfgs::lbfgs_parameter_t lbfgs_params;
    lbfgs_params.mem_size = 256;
    lbfgs_params.past = 10;
    lbfgs_params.g_epsilon = 0.0;
    lbfgs_params.min_step = 1e-32;
    lbfgs_params.delta = 1.0e-6;
    double minObjective;
    std::clock_t start = std::clock();
    auto ret = lbfgs::lbfgs_optimize(
        x_, minObjective, GlobTrajOpt::objectiveFunc, nullptr,
        nullptr, this, lbfgs_params);
    std::clock_t end = std::clock();

    std::cout << "\033[32m"
              << "ret: " << ret << ", minObjective: " << minObjective << ", iernum: " << iter_num << ", time: " << double(end - start) / CLOCKS_PER_SEC << "s" << "\033[0m"
              << std::endl;
    return ret;
  }

  inline void GlobTrajOpt::roundingState(Eigen::MatrixXd &state)
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

  inline void GlobTrajOpt::setBoundConds(const Eigen::MatrixXd &iniState,
                                         const Eigen::MatrixXd &finState)
  {
    jerkOpt_.reset(iniState, finState, N_);
  }

  double GlobTrajOpt::calTimeIntPenalty()
  {
    Eigen::Vector3d pos, vel, acc, jrk, snp;
    Eigen::Vector3d grad_tmp;
    double cost_tmp;
    static Eigen::Matrix<double, 6, 1> beta0, beta1, beta2, beta3, beta4;
    double s1, s2, s3, s4, s5;
    double step_vaj, omega, alpha;
    Eigen::Matrix<double, 6, 3> grad_Jv_to_c, grad_Ja_to_c, grad_Jj_to_c, grad_Jminthr_to_c, grad_Jmaxthr_to_c;
    Eigen::Matrix<double, 6, 3> integral_grad_Jv_to_c, integral_grad_Ja_to_c, integral_grad_Jj_to_c, integral_grad_Jminthr_to_c, integral_grad_Jmaxthr_to_c;
    double grad_Jv_to_T, grad_Ja_to_T, grad_Jj_to_T, grad_Jminthr_to_T, grad_Jmaxthr_to_T;
    double integral_grad_Jv_to_T, integral_grad_Ja_to_T, integral_grad_Jj_to_T, integral_grad_Jminthr_to_T, integral_grad_Jmaxthr_to_T;
    double integral_cost_v, integral_cost_a, integral_cost_j, integral_cost_minthr, integral_cost_maxthr;

    vel_pnt_ = 0.0;
    acc_pnt_ = 0.0;
    jrk_pnt_ = 0.0;
    thrmax_pnt_ = 0.0;
    thrmin_pnt_ = 0.0;
    double accumulated_dur(0.0);
    for (int i = 0; i < N_; ++i)
    {
      const auto &c = jerkOpt_.b.block<6, 3>(i * 6, 0);
      integral_grad_Jv_to_c.setZero();
      integral_grad_Ja_to_c.setZero();
      integral_grad_Jj_to_c.setZero();
      integral_grad_Jminthr_to_c.setZero();
      integral_grad_Jmaxthr_to_c.setZero();

      integral_grad_Jv_to_T = 0.0;
      integral_grad_Ja_to_T = 0.0;
      integral_grad_Jj_to_T = 0.0;
      integral_grad_Jminthr_to_T = 0.0;
      integral_grad_Jmaxthr_to_T = 0.0;

      integral_cost_v = 0.0;
      integral_cost_a = 0.0;
      integral_cost_j = 0.0;
      integral_cost_minthr = 0.0;
      integral_cost_maxthr = 0.0;

      step_vaj = jerkOpt_.T1(i) / config_.K;
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
        if (highDerivativeCostGrad(vel, 1, grad_tmp, cost_tmp))
        {
          grad_Jv_to_c = beta1 * grad_tmp.transpose();
          grad_Jv_to_T = alpha * grad_tmp.dot(acc);
          integral_grad_Jv_to_c += omega * grad_Jv_to_c;
          integral_grad_Jv_to_T +=
              omega * (cost_tmp / config_.K + step_vaj * grad_Jv_to_T);
          integral_cost_v += omega * cost_tmp;
        }
        if (highDerivativeCostGrad(acc, 2, grad_tmp, cost_tmp))
        {
          grad_Ja_to_c = beta2 * grad_tmp.transpose();
          grad_Ja_to_T = alpha * grad_tmp.dot(jrk);
          integral_grad_Ja_to_c += omega * grad_Ja_to_c;
          integral_grad_Ja_to_T +=
              omega * (cost_tmp / config_.K + step_vaj * grad_Ja_to_T);
          integral_cost_a += omega * cost_tmp;
        }
        if (highDerivativeCostGrad(jrk, 3, grad_tmp, cost_tmp))
        {
          grad_Jj_to_c = beta3 * grad_tmp.transpose();
          grad_Jj_to_T = alpha * grad_tmp.dot(snp);
          integral_grad_Jj_to_c += omega * grad_Jj_to_c;
          integral_grad_Jj_to_T +=
              omega * (cost_tmp / config_.K + step_vaj * grad_Jj_to_T);
          integral_cost_j += omega * cost_tmp;
        }

        // hzc thrust penalty
        Eigen::Vector3d thr = acc + Eigen::Vector3d(0.0, 0.0, 9.81);
        Eigen::Vector3d dSqrMagThr = 2 * thr;
        double fthr = thr.norm();
        double sqrMagThr = fthr * fthr;
        double violaThrl = thrAccMinSqr - sqrMagThr;
        double violaThrh = sqrMagThr - thrAccMaxSqr;
        // std::cout <<"thrAccMinSqr: "<< thrAccMinSqr<<" thrAccMaxSqr: "<<thrAccMaxSqr<<std::endl;
        if (violaThrl > 0.0)
        {
          //   std::cout << "violaThrl: "<<violaThrl<<std::endl;
          double violaPena(violaThrl), violaPenaGrad(1.0);
          /*
           grad = 2 * derivative * violaPenaGrad;
        cost += violaPena;
          */
          positiveSmoothedL1(violaThrl, config_.smoothEps, violaPena, violaPenaGrad);
          grad_Jminthr_to_c = -violaPenaGrad * beta2 * dSqrMagThr.transpose();
          grad_Jminthr_to_T = -violaPenaGrad * alpha * dSqrMagThr.transpose() * jrk;
          integral_grad_Jminthr_to_c += omega * grad_Jminthr_to_c;
          integral_grad_Jminthr_to_T +=
              omega * (violaPena / config_.K + step_vaj * grad_Jminthr_to_T);
          integral_cost_minthr += omega * violaPena;
        }

        if (violaThrh > 0.0)
        {
          double violaPena(violaThrh), violaPenaGrad(1.0);
          positiveSmoothedL1(violaThrh, config_.smoothEps, violaPena, violaPenaGrad);
          // std::cout <<"violaPena: "<<violaPena<<" violaThrh: "<<violaThrh<<" sqrMagThr: "<<sqrMagThr<< " acc: "<<acc.transpose()<<std::endl;

          grad_Jmaxthr_to_c = violaPenaGrad * beta2 * dSqrMagThr.transpose();
          grad_Jmaxthr_to_T = violaPenaGrad * alpha * dSqrMagThr.transpose() * jrk;
          integral_grad_Jmaxthr_to_c += omega * grad_Jmaxthr_to_c;
          integral_grad_Jmaxthr_to_T +=
              omega * (violaPena / config_.K + step_vaj * grad_Jmaxthr_to_T);
          integral_cost_maxthr += omega * violaPena;
        }
      }
      vel_pnt_ += step_vaj * integral_cost_v;
      acc_pnt_ += step_vaj * integral_cost_a;
      jrk_pnt_ += step_vaj * integral_cost_j;
      thrmin_pnt_ += step_vaj * integral_cost_minthr;
      thrmax_pnt_ += step_vaj * integral_cost_maxthr;
      jerkOpt_.gdC.block<6, 3>(i * 6, 0) +=
          step_vaj * (integral_grad_Jv_to_c * config_.pnlV +
                      integral_grad_Ja_to_c * config_.pnlA +
                      integral_grad_Jj_to_c * config_.pnlJ + integral_grad_Jminthr_to_c * config_.pnlThr +
                      integral_grad_Jmaxthr_to_c * config_.pnlThr);
      jerkOpt_.gdT(i) += integral_grad_Jv_to_T * config_.pnlV +
                         integral_grad_Ja_to_T * config_.pnlA +
                         integral_grad_Jj_to_T * config_.pnlJ + integral_grad_Jminthr_to_T * config_.pnlThr +
                         integral_grad_Jmaxthr_to_T * config_.pnlThr;

      accumulated_dur += jerkOpt_.T1(i);
    }
    vel_pnt_ *= config_.pnlV;
    acc_pnt_ *= config_.pnlA;
    jrk_pnt_ *= config_.pnlJ;
    thrmin_pnt_ *= config_.pnlThr;
    thrmax_pnt_ *= config_.pnlThr;
    return (vel_pnt_ + acc_pnt_ + jrk_pnt_ + thrmin_pnt_ + thrmax_pnt_);
  }

  /*
  * For penalty of higher derivatives like vel\acc\jrk\...
  f = max(v^2 - v_max^2, 0)
  */
  bool GlobTrajOpt::highDerivativeCostGrad(const Eigen::Vector3d &derivative,
                                           const double &squared_limit,
                                           Eigen::Vector3d &grad, double &cost)
  {
    double squared_diff = derivative.squaredNorm() - squared_limit;
    if (squared_diff > 0)
    {
      double violaPena(squared_diff), violaPenaGrad(1.0);
      positiveSmoothedL1(squared_diff, config_.smoothEps, violaPena,
                         violaPenaGrad);
      grad = 2 * derivative * violaPenaGrad;
      cost = violaPena;
      return true;
    }
    return false;
  }

  // tpye 1 for vel, 2 for acc, 3 for jerk
  bool GlobTrajOpt::highDerivativeCostGrad(const Eigen::Vector3d &derivative,
                                           const int &type, Eigen::Vector3d &grad,
                                           double &cost)
  {
    bool ret(false);
    cost = 0.0;
    grad.setZero();
    if (1 == type)
    {
      // velocity
      double sq_v_diff = derivative.squaredNorm() - squared_vel_limit_;
      if (sq_v_diff > DBL_EPSILON)
      {
        double violaPena(sq_v_diff), violaPenaGrad(1.0);
        positiveSmoothedL1(sq_v_diff, config_.smoothEps, violaPena, violaPenaGrad);
        grad = 2 * derivative * violaPenaGrad;
        cost += violaPena;
        ret = true;
      }
    }
    else if (2 == type)
    {
      // acceleration
      double sq_a_diff = derivative.squaredNorm() - squared_acc_limit_;
      if (sq_a_diff > DBL_EPSILON)
      {
        double violaPena(sq_a_diff), violaPenaGrad(1.0);
        positiveSmoothedL1(sq_a_diff, config_.smoothEps, violaPena, violaPenaGrad);
        grad = 2 * derivative * violaPenaGrad;
        cost += violaPena;
        ret = true;
      }
    }
    else if (3 == type)
    {
      // jerk
      double sq_j_diff = derivative.squaredNorm() - squared_jrk_limit_;
      if (sq_j_diff > DBL_EPSILON)
      {
        double violaPena(sq_j_diff), violaPenaGrad(1.0);
        positiveSmoothedL1(sq_j_diff, config_.smoothEps, violaPena, violaPenaGrad);
        grad = 2 * derivative * violaPenaGrad;
        cost += violaPena;
        ret = true;
      }
    }
    return ret;
  }

  // SECTION object function
  double GlobTrajOpt::objectiveFunc(void *ptrObj, const Eigen::VectorXd &x, Eigen::VectorXd &grad)
  {
    GlobTrajOpt &obj = *(GlobTrajOpt *)ptrObj;
    Eigen::Map<const Eigen::VectorXd> t(x.data(), obj.dim_t_);
    Eigen::Map<Eigen::VectorXd> gradt(grad.data(), obj.dim_t_);
    Eigen::VectorXd T(obj.N_);
    forwardT(t, T);
    // hzc
    Eigen::Map<const Eigen::VectorXd> endvel(x.data() + obj.dim_t_, 3);
    Eigen::Map<const Eigen::VectorXd> endacc(x.data() + obj.dim_t_ + 3, 3);
    Eigen::Map<Eigen::VectorXd> gradendv(grad.data() + obj.dim_t_, 3);
    Eigen::Map<Eigen::VectorXd> gradenda(grad.data() + obj.dim_t_ + 3, 3);
    Eigen::Matrix3d tailstate = obj.jerkOpt_.tailPVA;
    tailstate.col(1) = endvel;
    tailstate.col(2) = endacc;
    obj.jerkOpt_.resetTailCon(tailstate);

    obj.jerkOpt_.generate(obj.waypoints_, T);
    double cost = obj.jerkOpt_.getTrajJerkCost();
    obj.jerkOpt_.calGrads_CT();

    double time_itl_cost = obj.calTimeIntPenalty();
    cost += time_itl_cost;

    obj.jerkOpt_.calGrads_PT();

    /* add cost of rhoT * Ts and its grads to Ts */
    obj.jerkOpt_.gdT.array() += obj.config_.rhoTg;
    cost += obj.config_.rhoTg * T.sum();
    //   std::cout<<"time_cost: "<<obj.config_.rhoT * T.sum()<<std::endl;
    addLayerTGrad(t, obj.jerkOpt_.gdT, gradt);
    obj.iter_num += 1;

    gradendv = obj.jerkOpt_.gdTail.col(1);
    gradenda = obj.jerkOpt_.gdTail.col(2);

    gradendv.setZero();
    gradenda.setZero();

    return cost;
  }

  inline int GlobTrajOpt::earlyExit(void *ptrObj, const double *x,
                                    const double *grad, const double fx,
                                    const double xnorm, const double gnorm,
                                    const double step, int n, int k, int ls)
  {
    GlobTrajOpt &obj = *(GlobTrajOpt *)ptrObj;
    if (obj.vel_pnt_ <= 0.01 && obj.acc_pnt_ <= 0.01 && obj.jrk_pnt_ <= 0.01)
    {
      return 1;
    };
  }

  inline void GlobTrajOpt::forwardT(const Eigen::Ref<const Eigen::VectorXd> &t,
                                    Eigen::Ref<Eigen::VectorXd> vecT)
  {
    int M = t.size();
    for (int i = 0; i < M; ++i)
    {
      vecT(i) = expC2(t(i));
    }
    return;
  }

  inline void GlobTrajOpt::forwardTwithBound_Varied(const Eigen::Ref<const Eigen::VectorXd> &t,
                                                    Eigen::Ref<Eigen::VectorXd> vecT)
  {
    int M = t.size();
    for (int i = 0; i < M; ++i)
    {
      if (i < M - 1)
      {
        vecT(i) = expC2(t(i)) + 0.2 / (M - 1);
      }
      else
      {
        vecT(i) = expC2(t(i));
      }
    }
    return;
  }

  inline void GlobTrajOpt::backwardT(
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

  inline void GlobTrajOpt::backwardTwithBound_Varied(
      const Eigen::Ref<const Eigen::VectorXd> &vecT,
      Eigen::Ref<Eigen::VectorXd> t)
  {
    int M = t.size();
    for (int i = 0; i < M; ++i)
    {
      if (i < M - 1)
      {
        t(i) = logC2(vecT(i) - 0.2 / (M - 1));
      }
      else
      {
        t(i) = logC2(vecT(i));
      }
    }
    return;
  }

  inline void GlobTrajOpt::addLayerTGrad(
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

  // this is a smoothed C2 exact penalty
  // inputs x >0, output exact penalty and penalty derivates
  inline void GlobTrajOpt::positiveSmoothedL1(const double &x, const double &seps,
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

  inline double GlobTrajOpt::d_quasi_exp(const double &x)
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

  double GlobTrajOpt::expC2(double t)
  {
    return t > 0.0 ? ((0.5 * t + 1.0) * t + 1.0)
                   : 1.0 / ((0.5 * t - 1.0) * t + 1.0);
    // return 0.5 + t * t;
  }

  double GlobTrajOpt::logC2(double T)
  {
    return T > 1.0 ? (sqrt(2.0 * T - 1.0) - 1.0) : (1.0 - sqrt(2.0 / T - 1.0));
    // return sqrt(std::max(T - 0.5, 0.0));
  }

} // namespace traj_opt
