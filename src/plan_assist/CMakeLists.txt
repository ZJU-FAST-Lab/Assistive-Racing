cmake_minimum_required(VERSION 2.8.3)
project(plan_assist)

set(CMAKE_CXX_FLAGS "-std=c++14 ${CMAKE_CXX_FLAGS} -O3 -Wall")

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  visualization_msgs
  uav_utils
  message_filters
)

find_package(Eigen3 REQUIRED)
find_package(PCL 1.7 REQUIRED)

catkin_package(
 INCLUDE_DIRS include
 LIBRARIES plan_assist
)

include_directories( 
    SYSTEM 
    include 
    ${catkin_INCLUDE_DIRS}
    ${Eigen3_INCLUDE_DIRS} 
    ${PCL_INCLUDE_DIRS}
)

link_directories(${PCL_LIBRARY_DIRS})

add_library( plan_assist 
    src/RC.cpp
    src/risk_monitor.cpp
    src/corridor_helper.cpp
    )

target_link_libraries( plan_assist
    ${catkin_LIBRARIES} 
    ${PCL_LIBRARIES}
    )  
