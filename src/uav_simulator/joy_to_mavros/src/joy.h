#pragma once
#include <sensor_msgs/Joy.h>
#include <mavros_msgs/RCIn.h>
#include <ros/ros.h>
#include <vector>
#include <string>
using namespace std;
class JoyInput
{
public:
    sensor_msgs::Joy joy_msg;
    int joy_data_recived_flag = 0;
    void feed(sensor_msgs::JoyConstPtr msg);
};

class Joy
{
public:
    ros::Publisher joyRC_pub;
    int joy_channel_num;
    vector<int> joy_RC_joy_channel;
    vector<int> joy_RC_joy_reverse;
    JoyInput joyinput;
    double joy_freq_max;
    int calulateRC(int channel, int reverse);
    int remapRC(double input, int reverse);
    void pubRC();
    void config_from_ros_handle(const ros::NodeHandle &nh);

public:
private:
    template <typename TName, typename TVal>
    void read_essential_param(const ros::NodeHandle &nh, const TName &name, TVal &val)
    {
        if (nh.getParam(name, val))
        {
            // pass
        }
        else
        {
            ROS_ERROR_STREAM("Read param: " << name << " failed.");
            ROS_BREAK();
        }
    };
};