#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <ros/ros.h>

#include <string>
#include <vector>
#include <Eigen/Core>

namespace traj_opt
{
  struct Config
  {
    bool debug;
    int K;
    bool fix_end;
    double smoothEps;
    double rhoAccCmd; // zyh add
    double rhoVelCmd; // zyh add
    double rhoPosCmd; // zyh add
    double rhoErgBVP;
    double rhoT;
    double rhoTg;
    double rhoTdiff;
    double pnlTopo;
    double pnlP;
    double pnlV;
    double pnlA;
    double pnlJ;
    double pnlCone; // zyh add
    double dist0;   // zyh add
    double pnlThr;
    double pnlGate;
    double rhoClose;
    double pnlCap;
    double pnlSwarm;
    double pnlEqu;
    double rhoJCost;
    double rhoACost;
    double edt_dist_tolerance;
    double topo_dist_tolerance;
    double mutual_dist_tolerance;
    double tracking_dur;
    double capsule_radius;
    double capsule_up_t;
    double capsule_down_t;
    double capsule_dt;
    double Tau_cap;
    double Tau_obs;
    double Tau_swarm;
    double resolution;
    double max_vel;
    double max_acc;
    double max_jrk;
    double max_thracc;
    double min_thracc;
    double Tg_max;
    double Tg_min;
    double mass;
    double g;
    double hor_drag_coeff;
    double ver_drag_coeff;
    double par_drag_coeff;
    double speed_smooth_factor;
    double tilt_max;
    double body_rate_max;
    double thr_min;
    double thr_max;
    double rhoBodyRate;
    double rhoThr;
    double rhoTilt;
    double frontend_time_budget;
    int discretize_num_per_seg;
    double obs_radius;

    // zyh add
    double debug_sleep;
    double safe_traj_gen_time_;
    double pnlAligned;
    bool ifAligned;

    // yaw requirement
    double max_yaw_vel;
    double max_yaw_acc;
    double max_yaw_jrk;

    double pnlYawAligned;
    double rhoYawT;

    double rhoYawpos;
    double rhoYawvel;
    double rhoYawacc;

    double pnlyawV;
    double pnlyawA;
    double pnlyawJ;

    double rhoTup;

    double yawscalar;

    // progress requirement
    double T_max_vel;
    double T_max_acc;
    double T_max_jrk;

    bool dynamicRhoT;
    double rhoTProgress;
    double rhoTrajProgress;

    double pnltV;
    double pnltA;
    double pnltJ;

    double rhoTPositive;

    // fixed time integral pen
    double rc_t_lower_bound;
    double rc_K_fix;

    double yaw_t_lower_bound;
    double yaw_K_fix;

    // corridor search bound
    double search_forward_idx;
    double corridor_safe_dist;
    double corridor_edge_release_dist;

    // Load all parameters specified by ROS script
    void load(const ros::NodeHandle &nh)
    {
      nh.getParam("debug", debug);
      nh.getParam("K", K);
      nh.getParam("fix_end", fix_end);
      nh.getParam("smoothEps", smoothEps);

      nh.getParam("rhoAccCmd", rhoAccCmd);
      nh.getParam("rhoVelCmd", rhoVelCmd);
      nh.getParam("rhoPosCmd", rhoPosCmd);
      nh.getParam("pnlCone", pnlCone);
      nh.getParam("rhoErgBVP", rhoErgBVP);
      nh.getParam("rhoT", rhoT);
      nh.getParam("rhoTg", rhoTg);
      nh.getParam("rhoTdiff", rhoTdiff);
      nh.getParam("pnlTopo", pnlTopo);
      nh.getParam("pnlP", pnlP);
      nh.getParam("pnlV", pnlV);
      nh.getParam("pnlA", pnlA);
      nh.getParam("pnlJ", pnlJ);
      nh.getParam("pnlThr", pnlThr);
      nh.getParam("pnlGate", pnlGate);
      nh.getParam("rhoClose", rhoClose);
      nh.getParam("pnlCap", pnlCap);
      nh.getParam("pnlSwarm", pnlSwarm);
      nh.getParam("rhoJCost", rhoJCost);
      nh.getParam("rhoACost", rhoACost);
      nh.getParam("pnlEqu", pnlEqu);

      nh.getParam("edt_dist_tolerance", edt_dist_tolerance);
      nh.getParam("topo_dist_tolerance", topo_dist_tolerance);
      nh.getParam("mutual_dist_tolerance", mutual_dist_tolerance);

      nh.getParam("tracking_dur", tracking_dur);

      nh.getParam("capsule_radius", capsule_radius);
      nh.getParam("capsule_up_t", capsule_up_t);
      nh.getParam("capsule_down_t", capsule_down_t);
      nh.getParam("capsule_dt", capsule_dt);
      nh.getParam("Tau_cap", Tau_cap);
      nh.getParam("Tau_obs", Tau_obs);
      nh.getParam("Tau_swarm", Tau_swarm);

      nh.getParam("resolution", resolution);

      nh.getParam("max_vel", max_vel);
      nh.getParam("max_acc", max_acc);
      nh.getParam("max_jrk", max_jrk);
      nh.getParam("max_thracc", max_thracc);
      nh.getParam("min_thracc", min_thracc);

      nh.getParam("Tg_max", Tg_max);
      nh.getParam("Tg_min", Tg_min);

      nh.getParam("mass", mass);
      nh.getParam("g", g);
      nh.getParam("hor_drag_coeff", hor_drag_coeff);
      nh.getParam("ver_drag_coeff", ver_drag_coeff);
      nh.getParam("par_drag_coeff", par_drag_coeff);
      nh.getParam("speed_smooth_factor", speed_smooth_factor);
      nh.getParam("tilt_max", tilt_max);
      nh.getParam("body_rate_max", body_rate_max);
      nh.getParam("thr_min", thr_min);
      nh.getParam("thr_max", thr_max);
      nh.getParam("rhoBodyRate", rhoBodyRate);
      nh.getParam("rhoThr", rhoThr);
      nh.getParam("rhoTilt", rhoTilt);

      nh.getParam("frontend_time_budget", frontend_time_budget);
      nh.getParam("discretize_num_per_seg", discretize_num_per_seg);
      nh.getParam("obs_radius", obs_radius);

      // zyh add
      nh.getParam("dist0", dist0);
      nh.getParam("debug_sleep", debug_sleep);
      nh.getParam("safe_traj_gen_time", safe_traj_gen_time_);
      nh.getParam("pnlAligned", pnlAligned);
      nh.getParam("ifAligned", ifAligned);

      nh.getParam("max_yaw_vel", max_yaw_vel);
      nh.getParam("max_yaw_acc", max_yaw_acc);
      nh.getParam("max_yaw_jrk", max_yaw_jrk);

      nh.getParam("pnlYawAligned", pnlYawAligned);
      nh.getParam("rhoYawT", rhoYawT);

      nh.getParam("rhoYawpos", rhoYawpos);
      nh.getParam("rhoYawvel", rhoYawvel);
      nh.getParam("rhoYawacc", rhoYawacc);

      nh.getParam("pnlyawV", pnlyawV);
      nh.getParam("pnlyawA", pnlyawA);
      nh.getParam("pnlyawJ", pnlyawJ);

      nh.getParam("yawscalar", yawscalar);

      nh.getParam("T_max_vel", T_max_vel);
      nh.getParam("T_max_acc", T_max_acc);
      nh.getParam("T_max_jrk", T_max_jrk);

      nh.getParam("dynamicRhoT", dynamicRhoT);
      nh.getParam("rhoTProgress", rhoTProgress);
      nh.getParam("rhoTrajProgress", rhoTrajProgress);

      nh.getParam("pnltV", pnltV);
      nh.getParam("pnltA", pnltA);
      nh.getParam("pnltJ", pnltJ);

      nh.getParam("rhoTPositive", rhoTPositive);

      nh.getParam("rc_t_lower_bound", rc_t_lower_bound);
      nh.getParam("rc_K_fix", rc_K_fix);

      nh.getParam("yaw_t_lower_bound", yaw_t_lower_bound);
      nh.getParam("yaw_K_fix", yaw_K_fix);

      nh.getParam("search_forward_idx", search_forward_idx);
      nh.getParam("corridor_safe_dist", corridor_safe_dist);
      nh.getParam("corridor_edge_release_dist", corridor_edge_release_dist);

      nh.getParam("rhoTup", rhoTup);
    }
  };
} // namespace traj_opt

#endif