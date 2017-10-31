#include <string>
#include <iostream>
#include <ros/ros.h>
#include <ros/topic.h>
#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <sensor_msgs/Imu.h>

#include <boost/filesystem.hpp>

// Matlab header, this is needed to save mat files
// Note that we use the FindMatlab.cmake to get this
#include "mat.h"

using namespace std;

int main(int argc, char **argv) {

    // Debug message
    ROS_INFO("Starting up");

    // Check if there is a path to a dataset
    if(argc < 3) {
        ROS_ERROR("Error please specify a rosbag file");
        ROS_ERROR("Command Example: rosrun bagconvert bagconvert <rosbag> <topic>");
        return EXIT_FAILURE;
    }

    // Startup this node
    ros::init(argc, argv, "bagconvert");

    // Parse the input
    string pathBag = argv[1];
    string imuTopic = argv[2];

    // Get path
    boost::filesystem::path p(pathBag);
    string pathParent = p.parent_path().string();
    string pathMat;
    if(!pathParent.empty()) {
        pathMat = pathParent+"/"+p.stem().string()+".mat";
    } else {
        pathMat = p.stem().string()+".mat";
    }


    // Load rosbag here, and find messages we can play
    rosbag::Bag bag;
    bag.open(pathBag, rosbag::bagmode::Read);


    // We should load the bag as a view
    // Here we go from beginning of the bag to the end of the bag
    rosbag::View view(bag);

    // Debug
    ROS_INFO("BAG Path is: %s", pathBag.c_str());
    ROS_INFO("MAT Path is: %s", pathMat.c_str());
    ROS_INFO("Reading in rosbag file...");


    // Create the matlab mat file
    MATFile *pmat = matOpen(pathMat.c_str(), "w");
    if (pmat == NULL) {
        ROS_ERROR("Error could not create the mat file");
        return(EXIT_FAILURE);
    }

    // Our data vector
    vector<double> dataIMU = vector<double>();

    // Step through the rosbag and send to algo methods
    for (const rosbag::MessageInstance& m : view) {

        // Handle IMU message
        sensor_msgs::Imu::ConstPtr s1 = m.instantiate<sensor_msgs::Imu>();
        if (s1 != NULL && m.getTopic() == imuTopic) {
            dataIMU.push_back(m.getTime().toSec());
            dataIMU.push_back(s1->linear_acceleration.x);
            dataIMU.push_back(s1->linear_acceleration.y);
            dataIMU.push_back(s1->linear_acceleration.z);
            dataIMU.push_back(s1->angular_velocity.x);
            dataIMU.push_back(s1->angular_velocity.y);
            dataIMU.push_back(s1->angular_velocity.z);
        }

    }

    // Debug message
    ROS_INFO("Done processing bag");

    // ====================================================================
    // ==========              IMU DATA                  ==================
    // ====================================================================
    mxArray *pa1 = mxCreateDoubleMatrix(dataIMU.size()/7,7,mxREAL);
    if (pa1 == NULL) {
        printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
        printf("Unable to create mxArray.\n");
        return(EXIT_FAILURE);
    }
    // Correctly copy data over (column-wise)
    double* pt1 = mxGetPr(pa1);
    for(size_t i=0; i<dataIMU.size(); i+=7) {
        pt1[i/7] = dataIMU.at(i);
        pt1[(i + dataIMU.size())/7] = dataIMU.at(i+1);
        pt1[(i + 2*dataIMU.size())/7] = dataIMU.at(i+2);
        pt1[(i + 3*dataIMU.size())/7] = dataIMU.at(i+3);
        pt1[(i + 4*dataIMU.size())/7] = dataIMU.at(i+4);
        pt1[(i + 5*dataIMU.size())/7] = dataIMU.at(i+5);
        pt1[(i + 6*dataIMU.size())/7] = dataIMU.at(i+6);
    }
    // Add it to the matlab mat file
    int status = matPutVariable(pmat, "data_imu", pa1);
    if(status != 0) {
        printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
        return(EXIT_FAILURE);
    }
    // Cleanup
    mxDestroyArray(pa1);
    ROS_INFO("Done processing IMU data");

    // Close the mat file
    if (matClose(pmat) != 0) {
        ROS_ERROR("Error closing the mat file");
        return(EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
}


