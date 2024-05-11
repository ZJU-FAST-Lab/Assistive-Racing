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
#include <iostream>

#include "config.hpp"
#include "minco/minco_time_uniform.hpp"
#include "visualization/visualization.hpp"
// #define DEBUG

namespace traj_opt {

class UniformMincoStop {
   public:
    ros::NodeHandle nh_;
    Config config_;

    int waypoints_dim_, state_order_;
    int N_1_, dim_p_1_, dim_t_, dim_end_state_;
    Eigen::MatrixXd waypoints_1_, waypoints_2_;
    Eigen::VectorXd fixed_;

    int cost_function_cb_times_ = 0;

    Eigen::Vector3d gate_pos_;

    Eigen::MatrixX3d gdC_1_, gdC_2_;
    Eigen::Matrix3Xd gdP_1_, gdP_2_;
    Eigen::VectorXd gdT_, propagated_gdT_;

    // collision avoiding and dynamics paramters
    double squared_vel_limit_, squared_acc_limit_, squared_jrk_limit_;

    // Minimum Jerk Optimizer
    minco::uniform_T::UniformMinJerk jerkOpt_1_;
    // waypoints
    Eigen::VectorXd p_;
    // duration of each piece of the trajectory
    Eigen::VectorXd t_;

    // visualizer
    std::shared_ptr<visualization::Visualization> visPtr_;

    Eigen::VectorXd x_;
    Eigen::MatrixXd ini_state_, end_state_;

    // for penaly cost storage
    double vel_pnt_, acc_pnt_, jrk_pnt_, Cone_pnt_;

    Eigen::MatrixXd Safe_box_;

   public:
    UniformMincoStop(ros::NodeHandle& nh, Config& conf);
    ~UniformMincoStop() {}

    bool generate_traj(const Eigen::MatrixXd& iniState,
                       const Eigen::MatrixXd& finState,
                       const Eigen::MatrixXd& waypoints_1,
                       const double &dt_1,
                       Trajectory& traj);

    void setVisualizer(const std::shared_ptr<visualization::Visualization>& visPtr) {
        visPtr_ = visPtr;
    };

    void setSafeBox(const Eigen::MatrixXd safe_box){
        Safe_box_ = safe_box;
    }

   private:
    void setBoundConds(const Eigen::MatrixXd& iniState, const Eigen::MatrixXd& finState);
    void roundingState(Eigen::MatrixXd& state);
    int optimize(const double& delta = 1e-5);

    double calTimeIntPenalty(const int& N,
                             Eigen::MatrixX3d& gdC,
                             double& gdT,
                             double total_T,
                             const Eigen::MatrixX3d& coeffs);
    double calGatePenalty(const Eigen::MatrixXd& state,
                          Eigen::Ref<Eigen::MatrixXd> grad);

    bool safeHpolyCostGrad(const Eigen::Vector3d &pos,
                           Eigen::Vector3d &grad,
                           double &cost,
                           Eigen::MatrixXd hpolys);

    bool highDerivativeCostGrad_vel(const Eigen::Vector3d &derivative,
                                    Eigen::Vector3d &grad,
                                    double &cost);
    bool highDerivativeCostGrad_acc(const Eigen::Vector3d &derivative,
                                    Eigen::Vector3d &grad,
                                    double &cost);
    bool highDerivativeCostGrad_jrk(const Eigen::Vector3d &derivative,
                                    Eigen::Vector3d &grad,
                                    double &cost);
    static double objectiveFunc(void *ptrObj, const Eigen::VectorXd &x, Eigen::VectorXd& grad);

    static int Process(void* ptrObj,
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
    static void backwardT(const Eigen::Ref<const Eigen::VectorXd>& vecT, Eigen::Ref<Eigen::VectorXd> t);
    static void addLayerTGrad(const Eigen::Ref<const Eigen::VectorXd>& t,
                              const Eigen::Ref<const Eigen::VectorXd>& gradT,
                              Eigen::Ref<Eigen::VectorXd> gradt);
    
    static void addLayerPGrad(const Eigen::Ref<const Eigen::MatrixXd> &gradInPs,
                              Eigen::Ref<Eigen::MatrixXd> grad);

    // this is a smoothed C2 exact penalty
    // inputs x >0, output exact penalty and penalty derivates
    static void positiveSmoothedL1(const double& x,
                                   const double& seps,
                                   double& f,
                                   double& df);
    static double d_quasi_exp(const double& x);
    static double expC2(double t);
    static double logC2(double T);
};

}  // namespace traj_opt