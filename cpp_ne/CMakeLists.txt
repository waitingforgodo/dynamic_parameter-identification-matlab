cmake_minimum_required(VERSION 3.14)
project(robot_dynamic_torque_ne LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(robot_dyn_ne
    src/standard_regressor_marvin.cpp
    src/friction_regressor.cpp
    src/compute_torque.cpp
)
target_include_directories(robot_dyn_ne PUBLIC include)
