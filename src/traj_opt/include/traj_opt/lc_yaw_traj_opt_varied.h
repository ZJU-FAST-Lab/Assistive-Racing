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
// #include "minco/minco_time_uniform.hpp"
// #include "minco/minco_time_varied.hpp"
#include "minco/minco_yaw_time_varied.hpp"
#include "visualization/visualization.hpp"
#include "trajectory_utils/poly_traj_utils_7.hpp"

// #define DEBUG

namespace traj_opt {

class LocalYawOpt {
   public:
    ros::NodeHandle nh_;
    Config config_;

    int waypoints_dim_, state_order_;
    int N_, dim_p1_, dim_p2_,  dim_t1_, dim_t2_, dim_fin_;
    int N_1_, N_2_;
    Eigen::VectorXd waypoints1_, waypoints2_;
    Eigen::VectorXd yaw_waypoints1_, yaw_waypoints2_;
    Eigen::VectorXd fixed_;

    int cost_function_cb_times_ = 0;

    Eigen::Vector3d gate_pos_;

    Eigen::VectorXd gdC_1_, gdC_2_;
    Eigen::VectorXd gdP_1_, gdP_2_;
    Eigen::VectorXd gdT_1_, gdT_2_;

    Eigen::VectorXd gdT_, propagated_gdT_;

    // collision avoiding and dynamics paramters
    double squared_vel_limit_, squared_acc_limit_, squared_jrk_limit_;

    // Minimum Jerk Optimizer
    minco::MINCO_S31D YawOpt_;

    minco::MINCO_S31D lc_YawOpt_1_;
    minco::MINCO_S31D lc_YawOpt_2_;

    // waypoints
    Eigen::VectorXd p_;
    // duration of each piece of the trajectory
    Eigen::VectorXd t_;

    // visualizer
    std::shared_ptr<visualization::Visualization> visPtr_;

    Eigen::VectorXd x_;
    Eigen::VectorXd ini_state_, end_state_, mid_state_;

    // for penaly cost storage
    double vel_pnt_, acc_pnt_, jrk_pnt_, Cone_pnt_, pos_pnt_;

    //Assistance
    Trajectory1D Prim;
    Trajectory7 On_traj_;

    double rhoCmdAcc_, rhoCmdVel_, rhoCmdPos_;

    Eigen::MatrixXd Chpolys;

    int gate_pos_index;

    int K_fix_;
    double step_fix_;
    double t_lower_bound_;

    std::vector<Eigen::Vector3d> gate_pos_set_, gate_oritation_set_;
    std::vector<int> ids_;

    double Tsum_;
    double sumD1_, sumD2_;

    double forward_t_;
    
    double global_exp_scalar_;

    Eigen::Vector3d yaw_target_p;

   public:
    LocalYawOpt(ros::NodeHandle& nh, Config& conf);
    ~LocalYawOpt() {}

    bool generate_traj(const Eigen::VectorXd &iniState,
                       const Eigen::VectorXd &finState,
                       const Eigen::VectorXd &wps1,
                       const Eigen::VectorXd &D1,
                       Trajectory1D &traj);

    Trajectory1D get_yaw_minco_traj(const Eigen::VectorXd &iniState,
                                    const Eigen::VectorXd &finState,
                                    const Eigen::VectorXd &waypoints,
                                    const Eigen::VectorXd &dts);

    void setVisualizer(const std::shared_ptr<visualization::Visualization>& visPtr) {
        visPtr_ = visPtr;
    };

    void setRhoCmd(double rhoCmdAcc, double rhoCmdVel, double rhoCmdPos){
        rhoCmdAcc_ = rhoCmdAcc;
        rhoCmdVel_ = rhoCmdVel;
        rhoCmdPos_ = rhoCmdPos;
    }

    void setGateIndex(int idx){
        gate_pos_index = idx;
    }

    void setYawTarget(const Eigen::Vector3d global_yaw_target_p){
        yaw_target_p = global_yaw_target_p;
    }

    void setGatesPosAndVec(const std::vector<Eigen::Vector3d> &gate_pos_set, const std::vector<Eigen::Vector3d> &gate_oritation_set){
        gate_pos_set_ = gate_pos_set;
        gate_oritation_set_ = gate_oritation_set;
    }

    void setOnTraj(const Trajectory7 traj){
        On_traj_ = traj;
    }

    void setAssistTraj(const Trajectory1D primitive){
        Prim = primitive;
    }

    void setForwardT(const double t){
        forward_t_ = t;
    }

    // void setTimeForPrim(double t){
    //     t_for_prim = t;
    // }    

    void setConeParam(Eigen::Vector3d gate_pos, double lead, Eigen::Vector3d vec, double ratio);

    void sethpolys(Eigen::MatrixXd hpolys);

   private:
    void roundingState(Eigen::VectorXd& state);
    int optimize(const double& delta = 1e-5);

    bool SetOnTrajGateIndex();
    
    double calTimeIntPenalty(const int& N,
                             Eigen::VectorXd& gdC,
                             Eigen::VectorXd& gdT,
                             const Eigen::VectorXd& T,
                             const Eigen::VectorXd& coeffs,
                             const double t_start,
                             const double scalar) ;
    
    double calFixedTimeIntPenalty(const int& N,
                                 Eigen::VectorXd& gdC,
                                 Eigen::VectorXd& gdT,
                                 const Eigen::VectorXd& T,
                                 const Eigen::VectorXd& coeffs);

    bool UserIntentionPenalty(const double &acc,
                              const double &des_acc,
                              double &grad,
                              double &cost);

    bool VisibilityCostGrad(const Eigen::Vector3d &on_pos,
                            const Eigen::Vector3d &gate,
                            const double &phi,
                            double &grad,
                            Eigen::Vector3d &grad_on_pos,
                            double &cost);

    bool highDerivativeCostGrad_vel(const double &derivative,
                                    double &grad,
                                    double &cost);
    bool highDerivativeCostGrad_acc(const double &derivative,
                                    double &grad,
                                    double &cost);
    bool highDerivativeCostGrad_jrk(const double &derivative,
                                    double &grad,
                                    double &cost);

    static double objectiveFunc(void *ptrObj, const Eigen::VectorXd &x, Eigen::VectorXd& grad);

    static int Process(void *ptrObj, const Eigen::VectorXd &x,
                        const Eigen::VectorXd &grad, const double fx,
                        const double step, int k, int ls);

    void forwardT(const Eigen::Ref<const Eigen::VectorXd> &t, Eigen::Ref<Eigen::VectorXd> vecT, const double Ts);
    void backwardT(const Eigen::Ref<const Eigen::VectorXd>& vecT, Eigen::Ref<Eigen::VectorXd> t);

    
    static void addLayerTGrad(const Eigen::Ref<const Eigen::VectorXd> &t,
                              const double &Ts,
                              const Eigen::Ref<const Eigen::VectorXd> &gradT,
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