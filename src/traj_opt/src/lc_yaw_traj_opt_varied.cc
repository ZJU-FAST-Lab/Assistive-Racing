#include <traj_opt/lc_yaw_traj_opt_varied.h>

#include <traj_opt/chalk.hpp>
#include <optimization_opt/lbfgs_raw.hpp>

namespace traj_opt
{
  // time uniform
  LocalYawOpt::LocalYawOpt(ros::NodeHandle &nh, Config &conf)
      : nh_(nh), config_(conf)
  {
    squared_vel_limit_ = config_.max_yaw_vel * config_.max_yaw_vel;
    squared_acc_limit_ = config_.max_yaw_acc * config_.max_yaw_acc;
    squared_jrk_limit_ = config_.max_yaw_jrk * config_.max_yaw_jrk;
  }

  bool LocalYawOpt::generate_traj(const Eigen::VectorXd &iniState,
                                  const Eigen::VectorXd &finState,
                                  const Eigen::VectorXd &wps1,
                                  const Eigen::VectorXd &D1,
                                  Trajectory1D &traj)
  {
    ini_state_ = iniState;
    end_state_ = finState;

    N_1_ = D1.size();
    dim_t1_ = N_1_ - 1;

    sumD1_ = D1.sum();

    waypoints1_ = wps1;

    Eigen::VectorXd T1, T2;
    T1.resize(N_1_);
    T1 = D1;

    lc_YawOpt_1_.reset(N_1_);

    roundingState(ini_state_);
    roundingState(end_state_);

    dim_p1_ = wps1.size();
    dim_fin_ = finState.size();

    x_.resize(dim_p1_ + dim_t1_ + dim_fin_);
    Eigen::Map<Eigen::VectorXd> P1(x_.data(), dim_p1_);
    Eigen::Map<Eigen::VectorXd> t1(x_.data() + dim_p1_, dim_t1_);
    Eigen::Map<Eigen::VectorXd> finstate(x_.data() + dim_p1_ + dim_t1_, dim_fin_);

    P1 = wps1;

    backwardT(T1, t1); // time uniform unconstrained
    finstate = end_state_;

    t_lower_bound_ = config_.yaw_t_lower_bound;
    K_fix_ = config_.yaw_K_fix;
    step_fix_ = t_lower_bound_ / K_fix_;

    // std::cout<<"start optimize"<<std::endl;

    if (!SetOnTrajGateIndex())
      return false;

    // std::cout<<"current ids: "<<std::endl;
    // for(auto id : ids_){
    //   std::cout<<"[ "<<id<<" "<<gate_pos_set_[id].transpose()<<" ]"<<std::endl;
    // }
    // std::cout<<std::endl;

    int opt_ret = optimize();

    forwardT(t1, T1, sumD1_); // time varied with sum constraint

    if (opt_ret < 0)
    {
      return false;
    }

    // Time varied
    lc_YawOpt_1_.generate(ini_state_, finstate, P1, T1);

    traj = lc_YawOpt_1_.getTraj();

    return true;
  }

  bool LocalYawOpt::SetOnTrajGateIndex()
  {
    int num_yaw = N_1_ * config_.K + 1;
    double dt = sumD1_ / num_yaw;
    double s = forward_t_;
    Eigen::Vector3d pos, pos_ref, vel_ref, pos_ball;
    int id_tmp;
    Eigen::Vector3d gate_pos;

    ids_.clear();
    for (int i = 0; i < num_yaw; i++)
    {
      id_tmp = -1;
      pos = On_traj_.getPos(s);
      for (int j = gate_pos_index; j < int(gate_pos_set_.size()); j++)
      {
        // int k = j >= gate_pos_set_.size() ? j - gate_pos_set_.size() : j;
        gate_pos = gate_pos_set_[j];
        if ((gate_pos - pos).dot(gate_oritation_set_[j]) > 1e-2 && (gate_pos - pos).norm() > 1e-1)
        {
          id_tmp = j;
          break;
        }
      }

      ids_.push_back(id_tmp);
      s += dt;
    }
    return true;
  }

  Trajectory1D LocalYawOpt::get_yaw_minco_traj(const Eigen::VectorXd &iniState,
                                               const Eigen::VectorXd &finState,
                                               const Eigen::VectorXd &waypoints,
                                               const Eigen::VectorXd &dts)
  {

    ini_state_ = iniState;
    end_state_ = finState;

    N_ = waypoints.size() + 1;
    YawOpt_.reset(N_);

    roundingState(ini_state_);
    roundingState(end_state_);

    YawOpt_.generate(ini_state_, end_state_, waypoints, dts);

    return YawOpt_.getTraj();
  }

  int LocalYawOpt::optimize(const double &delta)
  {
    // Setup for L-BFGS solver
    lbfgs::lbfgs_parameter_t lbfgs_params;
    lbfgs_params.mem_size = 32;
    lbfgs_params.past = 3;
    lbfgs_params.g_epsilon = 0.0;
    lbfgs_params.min_step = 1e-64;
    lbfgs_params.delta = 1.0e-4;
    double minObjective;
    auto ret = lbfgs::lbfgs_optimize(
        x_, minObjective, LocalYawOpt::objectiveFunc, nullptr,
        //    LocalYawOpt::Process, this, lbfgs_params);//LocalYawOpt::Process
        nullptr, this, lbfgs_params);

    // std::cout << "\033[32m"
    //           << "[yaw_plan]ret: " << ret << ", minObjective: " << minObjective << "\033[0m"
    //           << std::endl;
    return ret;
  }

  inline void LocalYawOpt::roundingState(Eigen::VectorXd &state)
  {
    double tempNorm = abs(state(1));
    if (tempNorm > config_.max_yaw_vel)
    {
      state(1) *= (config_.max_yaw_vel - 0.001) / tempNorm;
      state(2) *= 0.0;
    }

    tempNorm = abs(state(2));
    if (tempNorm > config_.max_yaw_acc)
    {
      state(2) *= (config_.max_yaw_acc - 0.001) / tempNorm;
    }
  }

  double LocalYawOpt::calTimeIntPenalty(const int &N,
                                        Eigen::VectorXd &gdC,
                                        Eigen::VectorXd &gdT,
                                        const Eigen::VectorXd &T,
                                        const Eigen::VectorXd &coeffs,
                                        const double t_start,
                                        const double scalar)
  {
    //   std::cout<<"-------------------------"<<std::endl;
    double pos, vel, acc, jrk, snp;
    double grad_tmp;
    Eigen::Vector3d grad_on_pos;
    double cost_tmp;
    static Eigen::Matrix<double, 6, 1> beta0, beta1, beta2, beta3, beta4;
    double s1, s2, s3, s4, s5;
    double step_vaj, omega, alpha;
    Eigen::Matrix<double, 6, 1> grad_Jv_to_c, grad_Ja_to_c, grad_Jj_to_c, grad_Jp_to_c;
    Eigen::Matrix<double, 6, 1> integral_grad_Jv_to_c, integral_grad_Ja_to_c, integral_grad_Jj_to_c, integral_grad_Jp_to_c;
    double grad_Jv_to_T, grad_Ja_to_T, grad_Jj_to_T, grad_Jp_to_T;
    double integral_grad_Jv_to_T, integral_grad_Ja_to_T, integral_grad_Jj_to_T, integral_grad_Jp_to_T;
    double integral_cost_v, integral_cost_a, integral_cost_j, integral_cost_p;

    Eigen::Vector3d on_pos, on_vel;
    Eigen::Vector3d looking_gate;

    pos_pnt_ = 0.0;
    vel_pnt_ = 0.0;
    acc_pnt_ = 0.0;
    jrk_pnt_ = 0.0;
    Cone_pnt_ = 0.0;
    double on_t = t_start;
    int cnt = 0;
    int idx = 0;

    for (int i = 0; i < N; ++i)
    {
      const auto &c = coeffs.block<6, 1>(i * 6, 0);
      // const auto &on_c = On_traj_.getPiece(i).getCoeffMat();
      step_vaj = T[i] / config_.K;
      s1 = 0.0;

      integral_grad_Jv_to_c.setZero();
      integral_grad_Ja_to_c.setZero();
      integral_grad_Jj_to_c.setZero();
      integral_grad_Jp_to_c.setZero();

      integral_grad_Jv_to_T = 0.0;
      integral_grad_Ja_to_T = 0.0;
      integral_grad_Jj_to_T = 0.0;
      integral_grad_Jp_to_T = 0.0;

      integral_cost_v = 0.0;
      integral_cost_a = 0.0;
      integral_cost_j = 0.0;
      integral_cost_p = 0.0;

      for (int j = 0; j <= config_.K; ++j, s1 += step_vaj)
      {
        idx = i * config_.K + j;
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

        //   std::cout<<"on_t: "<<on_t<<std::endl;
        //   std::cout<<"pos: "<<pos<<"vel: "<<vel<<"acc: "<<acc<<"jrk: "<<jrk<<std::endl;
        on_pos = On_traj_.getPos(on_t);
        on_vel = On_traj_.getVel(on_t);

        alpha = (double)j / config_.K;
        omega = (j == 0 || j == config_.K) ? 0.5 : 1.0;
        if (highDerivativeCostGrad_vel(vel, grad_tmp, cost_tmp))
        {
          grad_Jv_to_c = beta1 * grad_tmp;
          grad_Jv_to_T = alpha * grad_tmp * acc;
          integral_grad_Jv_to_c += omega * grad_Jv_to_c;
          integral_grad_Jv_to_T +=
              omega * (cost_tmp / config_.K + step_vaj * grad_Jv_to_T);
          integral_cost_v += omega * cost_tmp;
        }
        if (highDerivativeCostGrad_acc(acc, grad_tmp, cost_tmp))
        {
          grad_Ja_to_c = beta2 * grad_tmp;
          grad_Ja_to_T = alpha * grad_tmp * jrk;
          integral_grad_Ja_to_c += omega * grad_Ja_to_c;
          integral_grad_Ja_to_T +=
              omega * (cost_tmp / config_.K + step_vaj * grad_Ja_to_T);
          integral_cost_a += omega * cost_tmp;
        }
        if (highDerivativeCostGrad_jrk(jrk, grad_tmp, cost_tmp))
        {
          grad_Jj_to_c = beta3 * grad_tmp;
          grad_Jj_to_T = alpha * grad_tmp * snp;
          integral_grad_Jj_to_c += omega * grad_Jj_to_c;
          integral_grad_Jj_to_T +=
              omega * (cost_tmp / config_.K + step_vaj * grad_Jj_to_T);
          integral_cost_j += omega * cost_tmp;
        }
        if (true)
        {
          // find which gate to look
          // for(int i = 0; i < gate_pos_set_.cols(); i++){
          //     gate_pos = gate_pos_set_.col(i);
          //     if((gate_pos - on_pos).dot(gate_oritation_set_.col(i)) > 1e-2 && (gate_pos - on_pos).norm() > 1e-1) break;
          // }
          looking_gate = ids_[idx] == -1 ? yaw_target_p : gate_pos_set_[ids_[idx]];
          if (VisibilityCostGrad(on_pos, looking_gate, pos, grad_tmp, grad_on_pos, cost_tmp))
          {
            // if(on_pos.norm() > 15){
            //         std::cout<<cnt<<": Trigger! gate_pos: "<<gate_pos_set_[ids_[idx]].transpose()<<" on_pos: "<<on_pos.transpose()<<std::endl;
            // }

            grad_Jp_to_c = beta0 * grad_tmp;
            grad_Jp_to_T = alpha * grad_tmp * vel + alpha * grad_on_pos.dot(on_vel);
            integral_grad_Jp_to_c += omega * grad_Jp_to_c;
            integral_grad_Jp_to_T +=
                omega * (cost_tmp / config_.K + step_vaj * grad_Jp_to_T);
            integral_cost_p += omega * cost_tmp;
          }
        }

        cnt++;
        if (j < config_.K)
          on_t += step_vaj;
      }

      vel_pnt_ += step_vaj * integral_cost_v;
      acc_pnt_ += step_vaj * integral_cost_a;
      jrk_pnt_ += step_vaj * integral_cost_j;
      pos_pnt_ += step_vaj * integral_cost_p;

      gdC.block<6, 1>(i * 6, 0) +=
          step_vaj * (integral_grad_Jv_to_c * config_.pnlyawV +
                      integral_grad_Ja_to_c * config_.pnlyawA +
                      integral_grad_Jj_to_c * config_.pnlyawJ +
                      integral_grad_Jp_to_c * config_.pnlYawAligned * scalar);

      gdT[i] += integral_grad_Jv_to_T * config_.pnlyawV +
                integral_grad_Ja_to_T * config_.pnlyawA +
                integral_grad_Jj_to_T * config_.pnlyawJ +
                integral_grad_Jp_to_T * config_.pnlYawAligned * scalar;
    }
    vel_pnt_ *= config_.pnlyawV;
    acc_pnt_ *= config_.pnlyawA;
    jrk_pnt_ *= config_.pnlyawJ;
    pos_pnt_ *= config_.pnlYawAligned * scalar;
    //   std::cout << "scalar : " << scalar << std::endl;

    // std::cout << "vel_pnt : " << vel_pnt_ << "\n"
    //           << "acc_pnt : " << acc_pnt_ << "\n"
    //           << "jrk_pnt : " << jrk_pnt_ << "\n"
    //           << "pos_pnt_ : " << pos_pnt_ << std::endl;

    return (vel_pnt_ + acc_pnt_ + jrk_pnt_ + Cone_pnt_ + pos_pnt_);
  }

  double LocalYawOpt::calFixedTimeIntPenalty(const int &N,
                                             Eigen::VectorXd &gdC,
                                             Eigen::VectorXd &gdT,
                                             const Eigen::VectorXd &T,
                                             const Eigen::VectorXd &coeffs)
  {
    double pos, vel, acc, jrk;
    double grad_tmp;
    double cost_tmp;
    static Eigen::Matrix<double, 6, 1> beta0, beta1, beta2, beta3, beta4;
    double s1, s2, s3, s4, s5;
    double omega;
    Eigen::Matrix<double, 6, 1> grad_Jp_to_c, grad_Jv_to_c, grad_Ja_to_c;
    Eigen::Matrix<double, 6, 1> integral_grad_Jp_to_c, integral_grad_Jv_to_c, integral_grad_Ja_to_c, integral_grad_Jfm_to_c;
    double grad_Jp_to_T, grad_Jv_to_T, grad_Ja_to_T;
    double integral_cost_p, integral_cost_v, integral_cost_a;

    int K_fix = K_fix_;
    double step_fix = step_fix_;

    int IdxPiece = 0;
    double t = 0.0;
    double t_pre = 0.0;

    integral_cost_p = 0;
    integral_cost_v = 0;
    integral_cost_a = 0;
    double rho = 1.0;

    for (int j = 0; j <= K_fix; ++j)
    {
      rho = exp2(-1.0 * double(j) / K_fix);
      // std::cout<<"j: "<<j<<std::endl;
      // std::cout<<"dT: "<<dT<<std::endl;
      // std::cout<<"t: "<<t<<std::endl;
      // std::cout<<"t_pre: "<<t_pre<<std::endl;
      while (t - t_pre > T(IdxPiece))
      {
        t_pre += T(IdxPiece);
        IdxPiece++;
      }

      // std::cout<<"IdxPiece: "<<IdxPiece<<std::endl;
      const auto &c = coeffs.block<6, 1>(IdxPiece * 6, 0);
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

      double acc_des = Prim.getAcc(t);
      double vel_des = Prim.getVel(t);
      double pos_des = Prim.getPos(t);
      // std::cout<<"t1: "<<t<<" acc_des: "<<acc_des.transpose()<<std::endl;
      omega = (j == 0 || j == K_fix) ? 0.5 : 1.0;

      if (UserIntentionPenalty(acc, acc_des, grad_tmp, cost_tmp))
      {
        // std::cout<<"tragger!acc: "<<acc.transpose()<<std::endl;
        grad_tmp = config_.rhoYawacc * grad_tmp;
        cost_tmp = config_.rhoYawacc * cost_tmp;
        grad_Ja_to_c = beta2 * grad_tmp;
        grad_Ja_to_T = grad_tmp * jrk;
        integral_grad_Ja_to_c = omega * step_fix * grad_Ja_to_c;
        if (IdxPiece > 0)
        {
          gdT.head(IdxPiece).array() += -rho * omega * step_fix * grad_Ja_to_T;
        }
        integral_cost_a += step_fix * omega * cost_tmp;
      }

      if (UserIntentionPenalty(vel, vel_des, grad_tmp, cost_tmp))
      {
        // std::cout<<"tragger!acc: "<<acc.transpose()<<std::endl;
        grad_tmp = config_.rhoYawvel * grad_tmp;
        cost_tmp = config_.rhoYawvel * cost_tmp;
        grad_Jv_to_c = beta1 * grad_tmp;
        grad_Jv_to_T = grad_tmp * acc;
        integral_grad_Jv_to_c = omega * step_fix * grad_Jv_to_c;
        if (IdxPiece > 0)
        {
          gdT.head(IdxPiece).array() += -rho * omega * step_fix * grad_Jv_to_T;
        }
        integral_cost_v += step_fix * omega * cost_tmp;
      }

      if (UserIntentionPenalty(pos, pos_des, grad_tmp, cost_tmp))
      {
        // std::cout<<"tragger!acc: "<<acc.transpose()<<std::endl;
        grad_tmp = config_.rhoYawpos * grad_tmp;
        cost_tmp = config_.rhoYawpos * cost_tmp;
        grad_Jp_to_c = beta0 * grad_tmp;
        grad_Jp_to_T = grad_tmp * vel;
        integral_grad_Jp_to_c = omega * step_fix * grad_Jp_to_c;
        if (IdxPiece > 0)
        {
          gdT.head(IdxPiece).array() += -rho * omega * step_fix * grad_Jp_to_T;
        }
        integral_cost_p += step_fix * omega * cost_tmp;
      }

      // std::cout<<"integral_grad_Ja_to_c: "<<integral_grad_Ja_to_c<<std::endl;
      // std::cout<<"integral_grad_Ja_to_T: "<<integral_grad_Ja_to_T<<std::endl;
      gdC.block<6, 1>(IdxPiece * 6, 0) += rho * (integral_grad_Ja_to_c + integral_grad_Jv_to_c + integral_grad_Jp_to_c);

      t += step_fix;
    }

    // std::cout<<"terminal successfully"<<std::endl;
    return integral_cost_a + integral_cost_v + integral_cost_p;
  }

  bool LocalYawOpt::UserIntentionPenalty(const double &val,
                                         const double &des_val,
                                         double &grad,
                                         double &cost)
  {
    bool ret(false);
    cost = 0.0;
    grad = 0.0;
    double pen = (val - des_val) * (val - des_val);
    if (pen > DBL_EPSILON)
    {
      double violaPena, violaPenaGrad;
      positiveSmoothedL1(pen, config_.smoothEps, violaPena, violaPenaGrad);
      grad = 2 * (val - des_val) * violaPenaGrad;
      cost += violaPena;
      ret = true;
    }

    return ret;
  }

  bool LocalYawOpt::VisibilityCostGrad(const Eigen::Vector3d &on_pos,
                                       const Eigen::Vector3d &gate,
                                       const double &phi,
                                       double &grad_phi,
                                       Eigen::Vector3d &grad_on_pos,
                                       double &cost)
  {
    // bool ret(false);
    // cost = 0.0;
    // grad = 0.0;
    // if((gate - on_pos).head(2).norm() < 0.1) return ret;
    // Eigen::Vector2d des_vec = (gate - on_pos).head(2);
    // double des_phi = atan2(des_vec[1], des_vec[0] + 1e-4);

    // double eval_val = (phi - des_phi) * (phi - des_phi);
    // if(eval_val > DBL_EPSILON){
    //     double violaPena, violaPenaGrad;
    //     positiveSmoothedL1(eval_val, config_.smoothEps, violaPena, violaPenaGrad);
    //     grad = 2 * (phi - des_phi) * violaPenaGrad;
    //     cost = violaPena;
    //     ret = true;
    //     // if(on_pos.norm() < 4.0){
    //     // std::cout<<"grad: "<<grad<<" ";
    //     // std::cout<<"cost: "<<cost<<" ";
    //     // std::cout<<"des_phi: "<<des_phi<<" ";
    //     // std::cout<<"phi: "<<phi<<std::endl;
    //     // }
    // }

    bool ret(false);
    cost = 0.0;
    grad_phi = 0.0;
    if ((gate - on_pos).head(2).norm() < 0.1)
      return ret;
    Eigen::Vector2d des_vec = (gate - on_pos).head(2).normalized();
    Eigen::Vector2d vec_phi(cos(phi), sin(phi));

    // double eval_val = 1.0 - vec_phi.dot(des_vec)/vec_phi.norm();
    double eval_val = 1.0 - vec_phi.dot(des_vec);
    if (eval_val > DBL_EPSILON)
    {
      double violaPena, violaPenaGrad;
      positiveSmoothedL1(eval_val, config_.smoothEps, violaPena, violaPenaGrad);
      // Eigen::Vector2d deval_to_vec_phi = - (1 /vec_phi.norm() * des_vec - vec_phi.dot(des_vec)/pow(vec_phi.norm(), 3) *  vec_phi);
      Eigen::Vector2d deval_to_vec_phi = -des_vec;
      Eigen::Vector2d dvec_phi_to_phi(-sin(phi), cos(phi));
      grad_phi = violaPenaGrad * deval_to_vec_phi.dot(dvec_phi_to_phi);

      Eigen::Vector3d des_vec_l = gate - on_pos;
      Eigen::Matrix3d Partialx_to_p = 1 / pow(des_vec_l.norm(), 3) * des_vec_l * des_vec_l.transpose() - 1 / des_vec_l.norm() * Eigen::Matrix3d::Identity();
      Eigen::MatrixXd Partialv_tp_x(3, 2);
      Partialv_tp_x << 1, 0,
          0, 1,
          0, 0;
      Eigen::Vector2d Partialeval_to_v = -vec_phi;
      grad_on_pos = Partialx_to_p * Partialv_tp_x * Partialeval_to_v * violaPenaGrad;
      // if(on_pos.norm() < 4.0){
      //     std::cout<<"-----------------------------------------------------------"<<std::endl;
      //     std::cout<<"phi: "<<phi<<std::endl;
      //     std::cout<<"des_vec: "<<des_vec.transpose()<<std::endl;
      //     std::cout<<"vec_phi: "<<vec_phi.transpose()<<std::endl;
      //     std::cout<<"violaPena: "<<violaPena<<std::endl;
      //     std::cout<<"violaPenaGrad: "<<violaPenaGrad<<std::endl;
      //     std::cout<<"deval_to_vec_phi: "<<deval_to_vec_phi.transpose()<<std::endl;
      //     std::cout<<"dvec_phi_to_phi: "<<dvec_phi_to_phi.transpose()<<std::endl;
      //     std::cout<<"grad: "<<grad<<std::endl;
      // }
      cost = violaPena;
      ret = true;
    }

    return ret;
  }

  /*
  * For penalty of higher derivatives like vel\acc\jrk\...
  f = max(v^2 - v_max^2, 0)
  */
  bool LocalYawOpt::highDerivativeCostGrad_vel(const double &derivative,
                                               double &grad,
                                               double &cost)
  {
    bool ret(false);
    cost = 0.0;
    grad = 0.0;
    double sq_v_diff = derivative * derivative - squared_vel_limit_;
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

  bool LocalYawOpt::highDerivativeCostGrad_acc(const double &derivative,
                                               double &grad,
                                               double &cost)
  {
    bool ret(false);
    cost = 0.0;
    grad = 0.0;
    // acceleration
    double sq_a_diff = derivative * derivative - squared_acc_limit_;
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

  bool LocalYawOpt::highDerivativeCostGrad_jrk(const double &derivative,
                                               double &grad,
                                               double &cost)
  {
    bool ret(false);
    cost = 0.0;
    grad = 0.0;
    // jerk
    double sq_j_diff = derivative * derivative - squared_jrk_limit_;
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

  double LocalYawOpt::objectiveFunc(void *ptrObj, const Eigen::VectorXd &x, Eigen::VectorXd &grad)
  {
    LocalYawOpt *obj = reinterpret_cast<LocalYawOpt *>(ptrObj);

    Eigen::Map<const Eigen::VectorXd> P1(x.data(), obj->dim_p1_);
    Eigen::Map<const Eigen::VectorXd> t1(x.data() + obj->dim_p1_, obj->dim_t1_);
    Eigen::Map<const Eigen::VectorXd> finstate(x.data() + obj->dim_p1_ + obj->dim_t1_, obj->dim_fin_);

    Eigen::Map<Eigen::VectorXd> gdP1(grad.data(), obj->dim_p1_);
    Eigen::Map<Eigen::VectorXd> gdt1(grad.data() + obj->dim_p1_, obj->dim_t1_);
    Eigen::Map<Eigen::VectorXd> gdfinstate(grad.data() + obj->dim_p1_ + obj->dim_t1_, obj->dim_fin_);

    // std::cout << "P1: " << P1.transpose() << std::endl;
    // std::cout << "t1: " << t1.transpose() << std::endl;
    // std::cout << "finstate: " << finstate.transpose() << std::endl;

    Eigen::VectorXd T1(obj->dim_t1_ + 1);
    obj->forwardT(t1, T1, obj->sumD1_);

    obj->lc_YawOpt_1_.generate(obj->ini_state_, finstate, P1, T1);
    // std::cout << "T1: " << T1.transpose() << std::endl;

    double cost = 0.0;
    cost += obj->lc_YawOpt_1_.getEnergy();

    // std::cout << "--------------------------------------------------" << std::endl;
    // std::cout << "energy_cost[lc_YawOpt_1_] : " << obj->lc_YawOpt_1_.getEnergy() << std::endl;
    // std::cout << "energy_cost[lc_YawOpt_2_] : " << obj->lc_YawOpt_2_.getEnergy() << std::endl;

    obj->lc_YawOpt_1_.calGrads_CT();

    // obj->lc_YawOpt_1_.gdC_Q.setZero();
    // obj->lc_YawOpt_2_.gdC_Q.setZero();
    // obj->lc_YawOpt_1_.gdT.setZero();
    // obj->lc_YawOpt_2_.gdT.setZero();

    obj->gdC_1_ = obj->lc_YawOpt_1_.gdC_Q;

    obj->gdT_1_ = obj->lc_YawOpt_1_.gdT;

    // std::cout << "calGrads_CT[gdC]: " << obj->YawOpt_.gdC_Q.transpose() << std::endl;
    // std::cout << "calGrads_CT[gdT]: " << obj->YawOpt_.gdT.transpose() << std::endl;

    // std::cout << "--------------------------------------------------" << std::endl;
    // std::cout << "start time integral penalty cal "<< std::endl;

    double t_sec = obj->forward_t_;
    // std::cout << "obj->forward_t_: "<<obj->forward_t_<< std::endl;
    obj->global_exp_scalar_ = 1.0;
    double time_itl_cost_1 = obj->calTimeIntPenalty(obj->N_1_, obj->gdC_1_, obj->gdT_1_, obj->lc_YawOpt_1_.T1, obj->lc_YawOpt_1_.b1, t_sec, 1.0);

    // std::cout << "TimeIntPenalty1 : " << "\n"
    //       << "vel_pnt : " << obj->vel_pnt_ << "\n"
    //       << "acc_pnt : " << obj->acc_pnt_ << "\n"
    //       << "jrk_pnt : " << obj->jrk_pnt_ << "\n"
    //       << "pos_pnt_ : " << obj->pos_pnt_ << std::endl;

    // t_sec = obj->forward_t_ + T1.sum();
    // double time_itl_cost_2 = obj->calTimeIntPenalty(obj->N_2_, obj->gdC_2_, obj->gdT_2_, obj->lc_YawOpt_2_.T1, obj->lc_YawOpt_2_.b1, obj->gate_pos_2_, t_sec, obj->config_.yawscalar);

    // std::cout << "TimeIntPenalty2 : " << "\n"
    //       << "vel_pnt : " << obj->vel_pnt_ << "\n"
    //       << "acc_pnt : " << obj->acc_pnt_ << "\n"
    //       << "jrk_pnt : " << obj->jrk_pnt_ << "\n"
    //       << "pos_pnt_ : " << obj->pos_pnt_ << std::endl;

    cost += time_itl_cost_1;

    // std::cout << "time_itl_cost_1 : " << time_itl_cost_1 << std::endl;

    // std::cout << "--------------------------------------------------" << std::endl;
    // std::cout << "start fixed time integral penalty cal "<< std::endl;

    double assist_cost = obj->calFixedTimeIntPenalty(obj->N_1_, obj->gdC_1_, obj->gdT_1_, obj->lc_YawOpt_1_.T1, obj->lc_YawOpt_1_.b1);
    cost += assist_cost;

    // std::cout << "assist_cost : " << assist_cost << std::endl;

    obj->lc_YawOpt_1_.gdC_Q = obj->gdC_1_;
    obj->lc_YawOpt_1_.gdT = obj->gdT_1_;

    obj->lc_YawOpt_1_.calGrads_PT();

    // std::cout << "calGrads_PT[Q]: " << obj->YawOpt_.gdQ.transpose() << std::endl;
    // std::cout << "calGrads_PT[T]: " << obj->YawOpt_.gdT.transpose() << std::endl;

    addLayerPGrad(obj->lc_YawOpt_1_.gdQ, gdP1);

    gdfinstate.setZero();
    gdfinstate += obj->lc_YawOpt_1_.gdTailQ;

    // std::cout << grad_mid_state.transpose() << std::endl;

    // cost += obj->config_.rhoT * (T1.sum() + T2.sum());
    // obj->lc_YawOpt_1_.gdT.array() += obj->config_.rhoYawT;
    // obj->lc_YawOpt_2_.gdT.array() += obj->config_.rhoYawT;

    addLayerTGrad(t1, obj->sumD1_, obj->lc_YawOpt_1_.gdT, gdt1);

    // std::cout << "gdP1: " << gdP1.transpose() << std::endl;
    // std::cout << "gdP2: " << gdP2.transpose() << std::endl;
    // std::cout << "gdt1: " << gdt1.transpose() << std::endl;
    // std::cout << "gdt2: " << gdt2.transpose() << std::endl;
    // std::cout << "gdmidstate: " << gdmidstate.transpose() << std::endl;

    // gdP1.setZero();
    // gdP2.setZero();
    gdt1.setZero();
    // gdt2.setZero();
    // gdmidstate.setZero();

    // grad_p_1.setZero();
    // grad_total_t.setZero();

    // grad_total_t.setZero();
    // grad_mid_state.col(0).setZero();

    return cost;
  }

  inline int LocalYawOpt::Process(void *ptrObj, const Eigen::VectorXd &x,
                                  const Eigen::VectorXd &grad, const double fx,
                                  const double step, int k, int ls)
  {
    LocalYawOpt *obj = reinterpret_cast<LocalYawOpt *>(ptrObj);

    if (obj->config_.debug == false)
    {
      return 0;
    }

    Eigen::Map<const Eigen::VectorXd> P1(x.data(), obj->dim_p1_);
    Eigen::Map<const Eigen::VectorXd> t1(x.data() + obj->dim_p1_, obj->dim_t1_);
    Eigen::Map<const Eigen::VectorXd> finstate(x.data() + obj->dim_p1_ + obj->dim_t1_, obj->dim_fin_);

    Eigen::VectorXd T1(obj->dim_t1_ + 1);
    obj->forwardT(t1, T1, obj->sumD1_);

    obj->lc_YawOpt_1_.generate(obj->ini_state_, finstate, P1, T1);

    Trajectory1D traj;
    traj = obj->lc_YawOpt_1_.getTraj();

    obj->visPtr_->visualize_1Dtraj(traj, "yaw_debug_traj");
    // Eigen::Map<const Eigen::VectorXd> tmp_grad(grad, obj->dim_p_1_ + obj->dim_p_2_ + obj->dim_t_ + obj->dim_mid_state_);
    // std::cout << "grad_norm = " << tmp_grad.norm() << std::endl;

    // Eigen::Map<const Eigen::MatrixXd> P1(x, obj->waypoints_dim_, obj->N_1_ - 1);
    // Eigen::Map<const Eigen::MatrixXd> P2(x + obj->dim_p_1_, obj->waypoints_dim_, obj->N_2_ - 1);
    // Eigen::Map<const Eigen::VectorXd> total_t(x + obj->dim_p_1_ + obj->dim_p_2_, obj->dim_t_);
    // Eigen::Map<const Eigen::MatrixXd> mid_state(x + obj->dim_p_1_ + obj->dim_p_2_ + obj->dim_t_, obj->waypoints_dim_, obj->state_order_);

    // Eigen::VectorXd total_T;
    // total_T.resize(obj->dim_t_);
    // // forwardT(total_t, total_T);
    // forwardTwithBound(total_t, total_T);

    // obj->jerkOpt_1_.resetTailCon(mid_state);
    // obj->jerkOpt_2_.resetHeadCon(mid_state);

    // obj->jerkOpt_1_.generate(P1, total_T[0]);
    // obj->jerkOpt_2_.generate(P2, total_T[1]);

    // Trajectory traj;
    // traj = obj->jerkOpt_1_.getTraj();
    // traj.append(obj->jerkOpt_2_.getTraj());

    // obj->visPtr_->visualize_traj(traj, "debug_traj");

    // double t = 0;
    // double dt = 0.01;
    // double t_scale = 3.0;
    // double duration = traj.getTotalDuration();

    // std::vector<Eigen::Vector3d> vn_list, vx_list, vy_list, vz_list;
    // std::vector<Eigen::Vector3d> an_list, ax_list, ay_list, az_list;
    // std::vector<Eigen::Vector3d> jn_list, jx_list, jy_list, jz_list;

    // for( t = 0; t < duration; t += dt)
    // {
    //     auto v = traj.getVel(t);
    //     auto a = traj.getAcc(t);
    //     auto j = traj.getJerk(t);

    //     vn_list.push_back(Eigen::Vector3d(t * t_scale, v.norm(), 0));
    //     vx_list.push_back(Eigen::Vector3d(t * t_scale, v[0], 0));
    //     vy_list.push_back(Eigen::Vector3d(t * t_scale, v[1], 0));
    //     vz_list.push_back(Eigen::Vector3d(t * t_scale, v[2], 0));
    //     an_list.push_back(Eigen::Vector3d(t * t_scale, a.norm(), 0));
    //     ax_list.push_back(Eigen::Vector3d(t * t_scale, a[0], 0));
    //     ay_list.push_back(Eigen::Vector3d(t * t_scale, a[1], 0));
    //     az_list.push_back(Eigen::Vector3d(t * t_scale, a[2], 0));
    //     jn_list.push_back(Eigen::Vector3d(t * t_scale, j.norm(), 0));
    //     jx_list.push_back(Eigen::Vector3d(t * t_scale, j[0], 0));
    //     jy_list.push_back(Eigen::Vector3d(t * t_scale, j[1], 0));
    //     jz_list.push_back(Eigen::Vector3d(t * t_scale, j[2], 0));
    // }

    // obj->visPtr_->visualize_path(vn_list, "vnorm");
    // obj->visPtr_->visualize_path(vx_list, "vx");
    // obj->visPtr_->visualize_path(vy_list, "vy");
    // obj->visPtr_->visualize_path(vz_list, "vz");
    // obj->visPtr_->visualize_path(an_list, "anorm");
    // obj->visPtr_->visualize_path(ax_list, "ax");
    // obj->visPtr_->visualize_path(ay_list, "ay");
    // obj->visPtr_->visualize_path(az_list, "az");
    // obj->visPtr_->visualize_path(jn_list, "jnorm");
    // obj->visPtr_->visualize_path(jx_list, "jx");
    // obj->visPtr_->visualize_path(jy_list, "jy");
    // obj->visPtr_->visualize_path(jz_list, "jz");

    ros::Duration(obj->config_.debug_sleep).sleep();

    return 0;
  }

  void LocalYawOpt::forwardT(const Eigen::Ref<const Eigen::VectorXd> &t,
                             Eigen::Ref<Eigen::VectorXd> vecT, const double Ts)
  {
    int M = t.size();
    for (int i = 0; i < M; ++i)
    {
      vecT(i) = expC2(t(i));
    }
    vecT(M) = 0.0;
    vecT /= 1.0 + vecT.sum();
    vecT(M) = 1.0 - vecT.sum();
    vecT *= Ts;
    return;
  }

  void LocalYawOpt::backwardT(const Eigen::Ref<const Eigen::VectorXd> &vecT,
                              Eigen::Ref<Eigen::VectorXd> t)
  {
    int M = t.size();
    Eigen::VectorXd vect(M);
    //   std::cout<<"[backwardT] M: "<<M<<std::endl;
    //   std::cout<<"[backwardT] vecT: "<<vecT.transpose()<<std::endl;
    double TM = vecT(M);
    for (int i = 0; i < M; i++)
    {
      vect[i] = vecT[i] / TM;
    }
    //   std::cout<<"[backwardT] vect: "<<vect.transpose()<<std::endl;
    //   t = vecT.head(M) / vecT(M);
    for (int i = 0; i < M; ++i)
    {
      t(i) = logC2(vect(i));
    }
    //   std::cout<<"[backwardT] t: "<<t.transpose()<<std::endl;
    return;
  }

  inline void LocalYawOpt::addLayerTGrad(const Eigen::Ref<const Eigen::VectorXd> &t,
                                         const double &Ts,
                                         const Eigen::Ref<const Eigen::VectorXd> &gradT,
                                         Eigen::Ref<Eigen::VectorXd> gradt)
  {
    int Ms1 = t.size();
    Eigen::VectorXd gFree = Ts * gradT.head(Ms1);
    double gTail = Ts * gradT(Ms1);
    Eigen::VectorXd dExpTau(Ms1);
    double expTauSum = 0.0, gFreeDotExpTau = 0.0;
    double denSqrt, expTau;
    for (int i = 0; i < Ms1; i++)
    {
      if (t(i) > 0)
      {
        expTau = (0.5 * t(i) + 1.0) * t(i) + 1.0;
        dExpTau(i) = t(i) + 1.0;
        expTauSum += expTau;
        gFreeDotExpTau += expTau * gFree(i);
      }
      else
      {
        denSqrt = (0.5 * t(i) - 1.0) * t(i) + 1.0;
        expTau = 1.0 / denSqrt;
        dExpTau(i) = (1.0 - t(i)) / (denSqrt * denSqrt);
        expTauSum += expTau;
        gFreeDotExpTau += expTau * gFree(i);
      }
    }
    denSqrt = expTauSum + 1.0;
    gradt = (gFree.array() - gTail) * dExpTau.array() / denSqrt -
            (gFreeDotExpTau - gTail * expTauSum) * dExpTau.array() / (denSqrt * denSqrt);
  }

  inline void LocalYawOpt::addLayerPGrad(
      const Eigen::Ref<const Eigen::MatrixXd> &gradInPs,
      Eigen::Ref<Eigen::MatrixXd> grad)
  {
    grad = gradInPs;
    return;
  }

  // this is a smoothed C2 exact penalty
  // inputs x >0, output exact penalty and penalty derivates
  inline void LocalYawOpt::positiveSmoothedL1(const double &x, const double &seps,
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

  inline double LocalYawOpt::d_quasi_exp(const double &x)
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

  double LocalYawOpt::expC2(double t)
  {
    return t > 0.0 ? ((0.5 * t + 1.0) * t + 1.0)
                   : 1.0 / ((0.5 * t - 1.0) * t + 1.0);
    // return 0.5 + t * t;
  }

  double LocalYawOpt::logC2(double T)
  {
    return T > 1.0 ? (sqrt(2.0 * T - 1.0) - 1.0) : (1.0 - sqrt(2.0 / T - 1.0));
    // return sqrt(std::max(T - 0.5, 0.0));
  }

  void LocalYawOpt::sethpolys(Eigen::MatrixXd hpolys)
  {
    Chpolys = hpolys;
  }

} // namespace traj_opt
