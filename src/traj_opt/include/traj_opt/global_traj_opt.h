/*
    MIT License

    Copyright (c) 2021 Hongkai Ye (kyle_yeh@163.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/
#pragma once

#include <ros/ros.h>

#include "config.hpp"
#include "minco/minco_time_varied.hpp"
#include "visualization/visualization.hpp"
// #define DEBUG

namespace traj_opt {

class GlobTrajOpt {
   public:
    ros::NodeHandle nh_;
    Config config_;

    bool fixed_end_;
    int N_, dim_t_;
    Eigen::MatrixXd waypoints_;

    // collision avoiding and dynamics paramters
    double squared_vel_limit_, squared_acc_limit_, squared_jrk_limit_,  thrAccMinSqr, thrAccMaxSqr;

    // Minimum Jerk Optimizer
    minco::varied_T::MinJerkOpt jerkOpt_;
    // waypoints
    Eigen::VectorXd p_;
    // duration of each piece of the trajectory
    Eigen::VectorXd t_;

    // visualizer
    std::shared_ptr<visualization::Visualization> visPtr_;

    // double* x_;
    Eigen::VectorXd x_;
    Eigen::MatrixXd ini_state_, end_state_;

    // for penaly cost storage
    double vel_pnt_, acc_pnt_, jrk_pnt_, thrmin_pnt_, thrmax_pnt_;

    double rhoCmdAcc_ = 1000000;

   public:
    GlobTrajOpt(ros::NodeHandle& nh, Config& conf);
    ~GlobTrajOpt() {}

    bool generate_traj(const Eigen::MatrixXd& iniState,
                       const Eigen::MatrixXd& finState,
                       const Eigen::MatrixXd& waypoints,
                       const Eigen::VectorXd& Ts,
                       Trajectory& traj);

    void setVisualizer(const std::shared_ptr<visualization::Visualization>& visPtr) {
        visPtr_ = visPtr;
    };

    void set_rhoT(double rhoT)
    {
        config_.rhoT = rhoT;
    }

   private:
    void setBoundConds(const Eigen::MatrixXd& iniState, const Eigen::MatrixXd& finState);
    void roundingState(Eigen::MatrixXd& state);
    int optimize(const double& delta = 1e-5);

    double calTimeIntPenalty();
    double calFixedTimeIntPenalty();

    bool highDerivativeCostGrad(const Eigen::Vector3d& derivative,
                                const double& squared_limit,
                                Eigen::Vector3d& grad,
                                double& cost);
    bool highDerivativeCostGrad(const Eigen::Vector3d& derivative,
                                const int& type,
                                Eigen::Vector3d& grad,
                                double& cost);
    bool UserIntentionPenalty(const Eigen::Vector3d &acc,
                              Eigen::Vector3d &des_acc,
                              Eigen::Vector3d &grad,
                              double &cost);

    static double objectiveFunc(void *ptrObj, const Eigen::VectorXd& x, Eigen::VectorXd &grad);
    static int earlyExit(void* ptrObj,
                         const double* x,
                         const double* grad,
                         const double fx,
                         const double xnorm,
                         const double gnorm,
                         const double step,
                         int n,
                         int k,
                         int ls);

    static void forwardT(const Eigen::Ref<const Eigen::VectorXd>& t, Eigen::Ref<Eigen::VectorXd> vecT);
    static void forwardTwithBound_Varied(const Eigen::Ref<const Eigen::VectorXd> &t, Eigen::Ref<Eigen::VectorXd> vecT);
    static void backwardT(const Eigen::Ref<const Eigen::VectorXd>& vecT, Eigen::Ref<Eigen::VectorXd> t);
    static void backwardTwithBound_Varied(const Eigen::Ref<const Eigen::VectorXd> &vecT, Eigen::Ref<Eigen::VectorXd> t);
    static void addLayerTGrad(const Eigen::Ref<const Eigen::VectorXd>& t,
                              const Eigen::Ref<const Eigen::VectorXd>& gradT,
                              Eigen::Ref<Eigen::VectorXd> gradt);

    // this is a smoothed C2 exact penalty
    // inputs x >0, output exact penalty and penalty derivates
    static void positiveSmoothedL1(const double& x,
                                   const double& seps,
                                   double& f,
                                   double& df);
    static double d_quasi_exp(const double& x);
    static double expC2(double t);
    static double logC2(double T);
    int iter_num;
};

}  // namespace traj_opt