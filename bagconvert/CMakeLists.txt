cmake_minimum_required(VERSION 2.8.3)

# Project name
project(bagconvert)

# Include our cmake files
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/)

# Include libraries
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem thread date_time)
find_package(Matlab 3 REQUIRED)


# Find catkin (the ROS build system)
find_package(catkin REQUIRED COMPONENTS roscpp nav_msgs std_msgs sensor_msgs nav_msgs rosbag)


# Describe catkin Project
catkin_package(
        DEPENDS Eigen3 Boost
        CATKIN_DEPENDS roscpp nav_msgs std_msgs sensor_msgs nav_msgs rosbag
        INCLUDE_DIRS src
)


# Try to compile with c++11
# http://stackoverflow.com/a/25836953
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# Enable compile optimizations
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fsee -fomit-frame-pointer -fno-signed-zeros -fno-math-errno -funroll-loops")

# Enable debug flags (use if you want to debug in gdb)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3  -Wall")



# Include our header files
include_directories(
        src
        ${EIGEN3_INCLUDE_DIR}
        ${Boost_INCLUDE_DIRS}
        ${MATLAB_INCLUDE_DIR}
        ${catkin_INCLUDE_DIRS}
)

# Set link libraries used by all binaries
list(APPEND parser_libraries
        ${Boost_LIBRARIES}
        ${MATLAB_LIBRARIES}
        ${MATLAB_MX_LIBRARY}
        ${MATLAB_MEX_LIBRARY}
        ${MATLAB_MAT_LIBRARY}
        ${MATLAB_ENG_LIBRARY}
        ${catkin_LIBRARIES}
)


##################################################
# Make binary for the offline reader
##################################################
add_executable(bagconvert
        src/main.cpp
)
target_link_libraries(bagconvert ${parser_libraries})

