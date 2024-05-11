#include <plan_manage/plan_manage.h>

int main(int argc, char **argv)
{
  ros::init(argc, argv, "planner_node");
  ros::NodeHandle nh("~");

  planner::PlanManager plan_manager(nh);

  ros::AsyncSpinner spinner(4);
  spinner.start();
  
  ros::waitForShutdown();

  return 0;
}


