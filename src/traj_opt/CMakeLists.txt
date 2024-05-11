cmake_minimum_required(VERSION 3.10.0)
project(traj_opt)

set(CMAKE_BUILD_TYPE "Release")

ADD_COMPILE_OPTIONS(-std=c++14 )

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")
# message( STATUS "OpenCV Include Dir: " ${OpenCV_INCLUDE_DIRS} )

find_package(Eigen3 REQUIRED)
find_package(PCL 1.7 REQUIRED)
find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
  visualization_msgs
  quadrotor_msgs
  uav_utils
  plan_assist
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
  DEPENDS Eigen 
  LIBRARIES traj_opt
)

add_library(traj_opt
    src/global_traj_opt.cc
    src/uniform_minco_progress.cc
    src/lc_yaw_traj_opt_varied.cc
    src/uniform_minco_emergency.cc
)

target_link_libraries(traj_opt
  ${catkin_LIBRARIES}
)

