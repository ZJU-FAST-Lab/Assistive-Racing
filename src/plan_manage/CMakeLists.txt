cmake_minimum_required(VERSION 3.10.0)
project(plan_manage)

set(CMAKE_BUILD_TYPE "Release")

ADD_COMPILE_OPTIONS(-std=c++14 )

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")
# message( STATUS "OpenCV Include Dir: " ${OpenCV_INCLUDE_DIRS} )

find_package(Eigen3 REQUIRED)
find_package(PCL 1.7 REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem)
find_package(catkin REQUIRED COMPONENTS
  roscpp
  roslib
  rospy
  std_msgs
  visualization_msgs
  quadrotor_msgs
  uav_utils
  tf2_ros
  decomp_ros_utils
  plan_assist
  traj_opt
)

include_directories(
    SYSTEM
    include
    ${catkin_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/include
    ${Eigen3_INCLUDE_DIRS}
    ${PCL_INCLUDE_DIRS}
)

catkin_package(
  INCLUDE_DIRS include
  DEPENDS Eigen Boost
)

add_executable(plan_manage
    src/main.cpp
    src/plan_manage.cpp
)

target_link_libraries(plan_manage
  ${catkin_LIBRARIES}
)

add_executable(fake_circle_node
    src/fake_circle_node.cpp
)

target_link_libraries(fake_circle_node
  ${catkin_LIBRARIES}
)

