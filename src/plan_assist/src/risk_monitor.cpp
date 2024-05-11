#include "risk_opt/risk_monitor.h"

static void positiveSmoothedL1(const double &x, const double &seps,
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

static double d_quasi_exp(const double &x)
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
static double expC2(double t)
{
  return t > 0.0 ? ((0.5 * t + 1.0) * t + 1.0)
                 : 1.0 / ((0.5 * t - 1.0) * t + 1.0);
}
static double logC2(double T)
{
  return T > 1.0 ? (sqrt(2.0 * T - 1.0) - 1.0) : (1.0 - sqrt(2.0 / T - 1.0));
}

namespace Risk_opt
{

  bool risk_monitor::ifAssist()
  {
    if (!if_ref_receive_)
      return false;
    if (!VelvecCheck(threshold_for_assist_))
      return true;
    return false;
  }

  bool risk_monitor::ifUserControl()
  {
    if (!if_ref_receive_)
      return false;
    if (VelvecCheck(threshold_for_prim_))
      return true;
    return false;
  }

  bool risk_monitor::SafetyCheck(Trajectory &emergency_traj)
  {
    if (!CorridorPtr_->ifCorridorReceived() || !if_state_receive_)
      return true;

    double init_T = 10;
    init_T = current_state_.col(1).norm() / risk_acc_max_;
    if (init_T < minSafeT_)
      init_T = minSafeT_ + 0.1;
    x_[0] = logC2(init_T - minSafeT_);

    bool ret = true;

    if (ret)
    {
      double T_fin = expC2(x_[0]) + minSafeT_;
      // std::cout<< "T_fin: "<<T_fin<<std::endl;
      Eigen::MatrixXd BoundMatrix(6, 6);
      BoundMatrix.setZero();
      BoundMatrix(0, 5) = 1.0;
      BoundMatrix(1, 4) = 1.0;
      BoundMatrix(2, 3) = 2.0;
      BoundMatrix(3, 4) = 1.0;
      BoundMatrix(3, 3) = 2 * T_fin;
      BoundMatrix(3, 2) = 3 * pow(T_fin, 2);
      BoundMatrix(3, 1) = 4 * pow(T_fin, 3);
      BoundMatrix(3, 0) = 5 * pow(T_fin, 4);
      BoundMatrix(4, 3) = 2.0;
      BoundMatrix(4, 2) = 6 * T_fin;
      BoundMatrix(4, 1) = 12 * pow(T_fin, 2);
      BoundMatrix(4, 0) = 20 * pow(T_fin, 3);
      BoundMatrix(5, 2) = 6;
      BoundMatrix(5, 1) = 24 * pow(T_fin, 1);
      BoundMatrix(5, 0) = 60 * pow(T_fin, 2);

      Eigen::MatrixXd BoundMatrix_inv = BoundMatrix.inverse();

      CoefficientMat c;
      c = (BoundMatrix_inv * b_).transpose();

      // std::cout<<"T_fin: "<<T_fin<<" c: "<<std::endl<<c.transpose()<<std::endl;

      emergency_traj = Trajectory(std::vector<Piece>{Piece(T_fin, c)});

      if (ifCollision(emergency_traj, c, T_fin))
        return false;

      // std::cout<<" No collision! "<<std::endl;
      return true;
    }
    return false;
  }

  bool risk_monitor::VelvecCheck(const double cos_threshold_)
  {
    Eigen::Vector3d pos = evaluate_traj.getPos(0.0);
    Eigen::Vector3d forward_direction = (next_gate_pos_ - pos).normalized();
    // std::cout<<"pos: "<<pos.transpose()<<std::endl;
    // std::cout<<"forward_direction: "<<forward_direction.transpose()<<std::endl;
    double forward_time_to_gate = evaluate_traj.getDurations().head(N_ + 1).sum();
    double current_t = (ros::Time::now() - eval_start_time_).toSec();
    for (double t = current_t; t < forward_time_to_gate; t += 0.01)
    {
      Eigen::Vector3d vel_tmp = evaluate_traj.getVel(t);
      if (vel_tmp.dot(forward_direction) / vel_tmp.norm() < cos_threshold_)
        return false;
    }
    return true;
  }

  double risk_monitor::SeekMinVelvec()
  {
    double cos_min = 1.0;
    Eigen::Vector3d pos = evaluate_traj.getPos(0.0);
    Eigen::Vector3d forward_direction = (next_gate_pos_ - pos).normalized();
    // std::cout<<"pos: "<<pos.transpose()<<std::endl;
    // std::cout<<"forward_direction: "<<forward_direction.transpose()<<std::endl;
    double forward_time_to_gate = evaluate_traj.getDurations().head(N_ + 1).sum();
    double current_t = (ros::Time::now() - eval_start_time_).toSec();
    for (double t = current_t; t < forward_time_to_gate; t += 0.01)
    {
      Eigen::Vector3d vel_tmp = evaluate_traj.getVel(t);
      if (vel_tmp.dot(forward_direction) / vel_tmp.norm() < cos_min)
      {
        cos_min = vel_tmp.dot(forward_direction) / vel_tmp.norm();
      }
    }
    return cos_min;
  }

  double risk_monitor::getDeviateCost(const ros::Time eval_t)
  {
    double t_for_eval = (eval_t - eval_start_time_).toSec();
    // int index = evaluate_traj.locatePieceIdx(t_for_eval);
    double eval_rest_time_ = evaluate_traj.getTotalDuration() - t_for_eval;

    double cum_cost = 0;
    double dt = 0.01;
    if (eval_rest_time_ > 1.0)
    {
      for (double t = 0.0; t < 1; t += dt)
      {
        cum_cost += pow(evaluate_traj.getAcc(t + t_for_eval).norm(), 2) * dt;
      }
    }
    else if (eval_rest_time_ > 0.0)
    {
      for (double t = 0.0; t < eval_rest_time_; t += dt)
      {
        cum_cost += pow(evaluate_traj.getAcc(t + t_for_eval).norm(), 2) * dt;
      }
    }
    acc_cost_ = cum_cost;

    return cum_cost;
  }

  double risk_monitor::time_evaluate(const ros::Time eval_t)
  {
    if (!if_ref_receive_ || !if_eval_receive_)
    {
      return -1;
    }
    double time_cost;
    double t_for_eval = (eval_t - eval_start_time_).toSec();

    Eigen::Vector3d pos_eval = evaluate_traj.getPos(t_for_eval);

    double eval_test_duration = evaluate_traj.getTotalDuration() - t_for_eval;
    double ref_test_duration = ref_traj.getTotalDuration() - getProgress(pos_eval);

    time_cost = (eval_test_duration - ref_test_duration) > 0 ? eval_test_duration - ref_test_duration : 0;

    time_cost_ = time_cost;

    return time_cost;
  }

  double risk_monitor::getProgress(const Eigen::Vector3d pos)
  {
    double eval_time = ref_traj.getDurations().head(N_ + 1).sum();
    double dist_min = INFINITY;
    Eigen::Vector3d pos_candidate;
    double time_candidate = 0.0;

    for (double t = 0; t < eval_time; t += 0.1)
    {
      Eigen::Vector3d pos_ref = ref_traj.getPos(t);
      double dist = (pos_ref - pos).norm();
      if (dist < dist_min)
      {
        dist_min = dist;
        pos_candidate = pos_ref;
        time_candidate = t;
      }
    }

    return time_candidate;
  }

  double risk_monitor::objectiveFunc(void *ptrObj, const Eigen::VectorXd &x, Eigen::VectorXd &grad)
  {
    risk_monitor *obj = reinterpret_cast<risk_monitor *>(ptrObj);

    const double t = x[0];
    // std::cout<<"t: "<<t<<std::endl;
    double gradT = 0;
    double cost = 0;
    double T = expC2(t) + obj->minSafeT_;

    Eigen::MatrixXd BoundMatrix(6, 6);
    BoundMatrix.setZero();
    BoundMatrix(0, 5) = 1.0;
    BoundMatrix(1, 4) = 1.0;
    BoundMatrix(2, 3) = 2.0;
    BoundMatrix(3, 4) = 1.0;
    BoundMatrix(3, 3) = 2 * T;
    BoundMatrix(3, 2) = 3 * pow(T, 2);
    BoundMatrix(3, 1) = 4 * pow(T, 3);
    BoundMatrix(3, 0) = 5 * pow(T, 4);
    BoundMatrix(4, 3) = 2.0;
    BoundMatrix(4, 2) = 6 * T;
    BoundMatrix(4, 1) = 12 * pow(T, 2);
    BoundMatrix(4, 0) = 20 * pow(T, 3);
    BoundMatrix(5, 2) = 6;
    BoundMatrix(5, 1) = 24 * pow(T, 1);
    BoundMatrix(5, 0) = 60 * pow(T, 2);

    Eigen::MatrixXd BoundMatrix_inv = BoundMatrix.inverse();

    CoefficientMat c;
    c = (BoundMatrix_inv * obj->b_).transpose();

    double dt = T / obj->K_;

    static Eigen::Matrix<double, 6, 1> beta0, beta1, beta2, beta3;
    double s1, s2, s3, s4, s5;
    Eigen::Vector3d pos, vel, acc, jrk;
    Eigen::Vector3d grad_tmp;
    double cost_tmp;
    double alpha, omega;

    s1 = 0;
    for (int i = 0; i <= obj->K_; i++)
    {
      s2 = s1 * s1;
      s3 = s2 * s1;
      s4 = s2 * s2;
      s5 = s4 * s1;
      beta0 << 1.0, s1, s2, s3, s4, s5;
      beta1 << 0.0, 1.0, 2.0 * s1, 3.0 * s2, 4.0 * s3, 5.0 * s4;
      beta2 << 0.0, 0.0, 2.0, 6.0 * s1, 12.0 * s2, 20.0 * s3;
      beta3 << 0.0, 0.0, 0.0, 6.0, 24.0 * s1, 60.0 * s2;
      beta0.reverse();
      beta1.reverse();
      beta2.reverse();
      beta3.reverse();

      pos = c * beta0;
      vel = c * beta1;
      acc = c * beta2;
      jrk = c * beta3;
      alpha = (double)i / obj->K_;
      omega = (i == 0 || i == obj->K_) ? 0.5 : 1.0;
      if (obj->MaxConstraint(vel, grad_tmp, cost_tmp, obj->risk_vel_max_))
      {
        gradT += obj->BigPen_ * dt * omega * alpha * grad_tmp.dot(acc);
        cost += obj->BigPen_ * dt * omega * cost_tmp;
      }
      if (obj->MaxConstraint(acc, grad_tmp, cost_tmp, obj->risk_acc_max_))
      {
        gradT += obj->BigPen_ * dt * omega * alpha * grad_tmp.dot(jrk);
        cost += obj->BigPen_ * dt * omega * cost_tmp;
      }

      s1 += dt;
    }

    cost += T;
    gradT += 1;

    grad[0] = gradT * d_quasi_exp(t);

    return cost;
  }

  bool risk_monitor::MaxConstraint(const Eigen::Vector3d &derivative,
                                   Eigen::Vector3d &grad,
                                   double &cost,
                                   const double &max_value)
  {
    bool ret(false);
    cost = 0.0;
    grad.setZero();
    double sq_v_diff = derivative.squaredNorm() - max_value * max_value;
    if (sq_v_diff > DBL_EPSILON)
    {
      double violaPena(sq_v_diff), violaPenaGrad(1.0);
      positiveSmoothedL1(sq_v_diff, 1.0e-2, violaPena, violaPenaGrad);
      grad = 2.0 * derivative * violaPenaGrad;
      cost += violaPena;
      ret = true;
    }

    return ret;
  }

  bool risk_monitor::ifCollision(const Trajectory &traj, const CoefficientMat &coeff, const double &duration)
  {
    // check the bound
    Eigen::VectorXd cx, cy, cz;
    Eigen::VectorXd check1;
    Eigen::Vector3d pos;
    int num = 10;
    double dt = duration / num;

    double s = 0.0;
    int id_last = current_id_;
    for (int j = 0; j <= num; j++)
    {
      Eigen::Vector3d pos_end = traj.getPos(s);
      int id1, id2;

      id1 = CorridorPtr_->SafeCheckHpolyForward(pos_end, id_last);
      id2 = CorridorPtr_->SafeCheckHpolyBackward(pos_end, id_last);
      if (!(id1 == -1) || !(id2 == -1))
      {
        if (id1 != -1)
        {
          id_last = id1;
          s += dt;

          continue;
        }
        else
        {
          id_last = id2;
          s += dt;
          continue;
        }
      }
      return true;
    }

    return false;
  }

};
