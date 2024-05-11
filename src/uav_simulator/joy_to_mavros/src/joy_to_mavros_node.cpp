#include "joy.h"
#include <sensor_msgs/Joy.h>
int main(int argc, char *argv[])
{
  ros::init(argc, argv, "joy_to_mavros");
  ros::NodeHandle nh("~");

  Joy joy;
  joy.config_from_ros_handle(nh);

  ros::Subscriber joy_sub =
      nh.subscribe<sensor_msgs::Joy>("/joy",
                                     100,
                                     boost::bind(&JoyInput::feed, &joy.joyinput, _1),
                                     ros::VoidConstPtr(),
                                     ros::TransportHints().tcpNoDelay());

  joy.joyRC_pub = nh.advertise<mavros_msgs::RCIn>("in",
                                                  100);
  ros::Rate r(joy.joy_freq_max);

  while (ros::ok())
  {
    r.sleep();
    ros::spinOnce();
    joy.pubRC();
  }
}
