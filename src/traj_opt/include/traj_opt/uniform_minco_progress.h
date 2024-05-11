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
#include <random>
#include <cmath>

#include "config.hpp"
#include "minco/minco_time_uniform.hpp"
#include "minco/minco_yaw_time_uniform.hpp"
#include "minco/minco_yaw_time_varied.hpp"
#include "visualization/visualization.hpp"
#include "corridor_helper/corridor_helper.h"
// #include "plan_env/sdf_map.h"
// #define DEBUG

namespace traj_opt {

struct PolyLine{
    Eigen::Matrix<double, 1, 6> line;
    std::vector<Eigen::Vector3d> points;  // intersecting point
    std::vector<int> intersecting_line_id;  // 
};

typedef Eigen::MatrixX4d Polyhedron;

class UniformMincoProgress {
   public:
    ros::NodeHandle nh_;
    Config config_;

    int waypoints_dim_, state_order_;
    int N_1_, dim_p_1_, dim_t_p_, dim_end_state_, dim_tend_state_;
    Eigen::MatrixXd waypoints_1_;
    Eigen::VectorXd Twaypoints_1_;

    int cost_function_cb_times_ = 0;

    Eigen::Vector3d gate_pos_;

    Eigen::MatrixX3d gdC_1_;
    Eigen::Matrix3Xd gdP_1_;
    Eigen::VectorXd gdC_t_;
    Eigen::VectorXd gdP_t_;
    // Eigen::VectorXd gdT_, propagated_gdT_;

    // collision avoiding and dynamics paramters
    double squared_vel_limit_, squared_acc_limit_, squared_jrk_limit_;

    // Minimum Jerk Optimizer
    // minco::uniform_T::UniformMinJerk jerkOpt_1_;
    minco::uniform_T::UniformMinSnap snapOpt_1_;
    minco::uniform_T::UniformMinJerk1D Timeopt_1_;
    minco::MINCO_S31D Timeopt_2_;

    // waypoints
    Eigen::VectorXd p_;
    // duration of each piece of the trajectory
    Eigen::VectorXd t_;
    // total duration
    double total_t_;

    // visualizer
    std::shared_ptr<visualization::Visualization> visPtr_;
    std::shared_ptr<corridor_opt::CorridorHelper> CorridorPtr_;

    Eigen::VectorXd x_;
    Eigen::MatrixXd ini_state_, mid_state_, end_state_;
    Eigen::VectorXd ini_T_, fin_T_;

    // for penaly cost storage
    double vel_pnt_, acc_pnt_, jrk_pnt_, Cone_pnt_;
    double t_vel_pnt_, t_acc_pnt_, t_jrk_pnt_;
    double t_progress_cost_;
    double pos_progress_cost_;
    double t_positive_pnt_;
    double t_bound_pnt_;

    double Rc_pos_pnt_, Rc_vel_pnt_, Rc_acc_pnt_;

    //Assistance
    Trajectory7 Prim;
    double dist0_;

    double resolution_inv_;

    //reference
    Trajectory ref_traj_;

    // fixed time integral penalty
    int K_fix_;
    double step_fix_;
    double t_lower_bound;

    std::vector<int> ids_;
    std::vector<int> d_set_;
    int current_id_;

    //debug for accelerating the computation
    double perperation_t1_, perperation_t2_, perperation_t3_, optimization_t_, iteration_;
    std::vector<Eigen::MatrixXd> Hpolys_in_;
    std::vector<std::pair<double, double>> time_searial_;

    Eigen::Vector3d vel_cmd_;
    double rhoTprogress_;

   public:
    UniformMincoProgress(ros::NodeHandle& nh, Config& conf);
    ~UniformMincoProgress() {}

    void setInitialState(const Eigen::MatrixXd &iniState, 
                         const Eigen::MatrixXd &finState,
                         const Eigen::VectorXd &initT,
                         const Eigen::VectorXd &finT,
                         const int N);

    bool generate_traj(const Eigen::MatrixXd &Poswaypoints,
                       const Eigen::VectorXd &Twaypoints,
                       const double &Total_T,
                       Trajectory7 &traj,
                       Trajectory1D &time_traj,
                       Trajectory7 &primitive);

    void setVisualizer(const std::shared_ptr<visualization::Visualization>& visPtr) {
        visPtr_ = visPtr;
    };

    void setCorridorHelper(const std::shared_ptr<corridor_opt::CorridorHelper>& CorridorPtr) {
        CorridorPtr_ = CorridorPtr;
    };    

    void setRefTrajectory(const Trajectory &traj){
        ref_traj_ = traj;
    }

    void setCurrentCorridorId(const int idx){
        current_id_ = idx;
    }

    void setTimeSerial(const std::vector<std::pair<double, double>> time_searial){
        time_searial_ = time_searial;
    }

    void setVelCmd(const Eigen::Vector3d vel_cmd){
        vel_cmd_  = vel_cmd;
    }

    double getDeviateCost(){
        return pos_progress_cost_/config_.rhoTrajProgress;
    }

    double getRCControlCost(){
        return Rc_pos_pnt_/config_.rhoPosCmd + Rc_vel_pnt_/config_.rhoVelCmd + Rc_acc_pnt_/config_.rhoAccCmd;
    }

    double getCorridorPenalty(){
        return Cone_pnt_;
    }

    double getRhoTProgress(){
        return rhoTprogress_;
    }

    double getIteration(){
        return iteration_;       
    }

    std::vector<Eigen::MatrixXd> getHpolysIn(){
        return Hpolys_in_;       
    }

    std::vector<double> getComputationTimes(){
        std::vector<double> T;
        T.push_back(perperation_t1_);
        T.push_back(perperation_t2_);
        T.push_back(perperation_t3_);
        T.push_back(optimization_t_);

        return T;
    }

   private:
    void setBoundConds(const Eigen::MatrixXd &iniState,
                       const Eigen::MatrixXd &finState,
                       const Eigen::VectorXd &iniT,
                       const Eigen::VectorXd &finT);

    void roundingState(Eigen::MatrixXd &state, Eigen::VectorXd &T);

    void calculatePrimSafeHpolys();

    Eigen::Vector3d findNeareastPos(const Eigen::Vector3d pos, const double heuristic_t, int id_last, double &t_best);

    Eigen::Vector3d findNeareastPos_numerical(const Eigen::Vector3d pos, const double heuristic_t, int id_last, double &t_best);

    void calculateSafeHpolys(const Trajectory7 pos_traj, const Trajectory1D time_traj);

    int optimize(const double& delta = 1e-5);

    double calTimeIntPenalty(const int& N,
                            Eigen::MatrixX3d& gdC,
                            Eigen::VectorXd& gdCt); 

    double calFixedTimeIntPenalty(const int& N,
                                  Eigen::MatrixX3d& gdC);

    bool UserIntentionPenalty(const Eigen::Vector3d &acc,
                              Eigen::Vector3d &des_acc,
                              Eigen::Vector3d &grad,
                              double &cost);       

    double calGatePenalty(const Eigen::MatrixXd& state,
                          Eigen::Ref<Eigen::MatrixXd> grad);

    bool ProgressCostGrad(const Eigen::Vector3d &pos,
                          const Eigen::Vector3d &ref_p,
                          Eigen::Vector3d &gradp,
                          Eigen::Vector3d &gradrefp,
                          double &cost);

    bool TimeboundCostGrad(const double &t,
                           double &grad,
                           double &cost);

    bool safeHpolyCostGrad(const Eigen::Vector3d &pos,
                           Eigen::Vector3d &grad,
                           double &cost,
                           Eigen::MatrixXd hpolys);

    bool NormalMaxConstriant(const Eigen::Vector3d &derivative,
                             Eigen::Vector3d &grad,
                             double &cost,
                             const double max_limit);

    bool NormalMaxConstriant1D(const double &derivative,
                               double &grad,
                               double &cost,
                               const double max_limit);

    bool PostiveTimeConstraint(const double &t,
                               double &grad,
                               double &cost); 

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
    
    static int earlyExit(void *ptrObj, const Eigen::VectorXd &x,
                         const Eigen::VectorXd &grad, const double fx,
                         const double step, int k, int ls);

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