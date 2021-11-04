/**
* This file is part of ORB-SLAM2.
*
* Copyright (C) 2014-2016 RaÃºl Mur-Artal <raulmur at unizar dot es> (University of Zaragoza)
* For more information see <https://github.com/raulmur/ORB_SLAM2>
*
* ORB-SLAM2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with ORB-SLAM2. If not, see <http://www.gnu.org/licenses/>.
*/


#include<iostream>
#include<algorithm>
#include<fstream>
#include<chrono>

#include<opencv2/core/core.hpp>

#include<System.h>

//#include"cuda_runtime.h"

using namespace std;

void LoadImages(const string &strAssociationFilename, string& file_type, string& basic_rgb, string& basic_depth, vector<string> &vstrImageFilenamesRGB,
                vector<string> &vstrImageFilenamesD, vector<double> &vTimestamps);

int main(int argc, char **argv)
{
    if(argc != 6)
    {
        cerr << endl << "Usage: ./rgbd_tsukuba path_to_vocabulary path_to_settings path_to_sequence path_to_association" << endl;
        return 1;
    }

    cudaFree(0);

    // Retrieve paths to images
    vector<string> vstrImageFilenamesRGB;
    vector<string> vstrImageFilenamesD;
    vector<double> vTimestamps;
    string file_type = string(argv[3]);
    string basic_rgb_add= string(argv[4]);
    string strAssociationFilename = string(argv[5]);
    strAssociationFilename = basic_rgb_add  + "/" + file_type + "/times.txt";
    basic_rgb_add += "/" + file_type + "/left";
    string basic_depth_add = string(argv[5]);
    LoadImages(strAssociationFilename, file_type, basic_rgb_add, basic_depth_add, vstrImageFilenamesRGB, vstrImageFilenamesD, vTimestamps);

    //std::cout<<"size of rgb, depth and timestamp "<<vstrImageFilenamesRGB.size()<<","<<vstrImageFilenamesD.size()<<","<<vTimestamps.size()<<std::endl;

    // Check consistency in the number of images and depthmaps
    int nImages = vstrImageFilenamesRGB.size();
    if(vstrImageFilenamesRGB.empty())
    {
        cerr << endl << "No images found in provided path." << endl;
        return 1;
    }
    else if(vstrImageFilenamesD.size()!=vstrImageFilenamesRGB.size())
    {
        cerr << endl << "Different number of images for rgb and depth." << endl;
        return 1;
    }

    // // Create SLAM system. It initializes all system threads and gets ready to process frames.
    ORB_SLAM2::System SLAM(argv[1],argv[2],ORB_SLAM2::System::RGBD,true);

    // // Vector for tracking time statistics
    vector<float> vTimesTrack;
    vTimesTrack.resize(nImages);

    cout << endl << "-------" << endl;
    cout << "Start processing sequence ..." << endl;
    cout << "Images in the sequence: " << nImages << endl << endl;

    // // Main loop
    cv::Mat imRGB, imD;
    cv::FileStorage fs;
    
    for(int ni=0; ni<nImages; ni++)
    {
        // Read image and depthmap from file
        imRGB = cv::imread(vstrImageFilenamesRGB[ni],CV_LOAD_IMAGE_UNCHANGED);
        fs.open(vstrImageFilenamesD[ni], cv::FileStorage::READ);
        fs["depth"]>>imD;

        double tframe = vTimestamps[ni];

        if(imRGB.empty())
        {
            cerr << endl << "Failed to load image at: "
                 <<vstrImageFilenamesRGB[ni] << endl;
            return 1;
        }
        
#ifdef COMPILEDWITHC11
        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
#else
        std::chrono::monotonic_clock::time_point t1 = std::chrono::monotonic_clock::now();
#endif

        // Pass the image to the SLAM system
        SLAM.TrackRGBD(imRGB,imD,tframe);

#ifdef COMPILEDWITHC11
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
#else
        std::chrono::monotonic_clock::time_point t2 = std::chrono::monotonic_clock::now();
#endif

        double ttrack= std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();

        vTimesTrack[ni]=ttrack;

        // Wait to load the next frame
        double T=0;
        if(ni<nImages-1)
            T = vTimestamps[ni+1]-tframe;
        else if(ni>0)
            T = tframe-vTimestamps[ni-1];

        if(ttrack<T)
            usleep((T-ttrack)*1e6);

    }

    // // Stop all threads
    SLAM.Shutdown();
    fs.release();

    // Tracking time statistics
    sort(vTimesTrack.begin(),vTimesTrack.end());
    float totaltime = 0;
    for(int ni=0; ni<nImages; ni++)
    {
        totaltime+=vTimesTrack[ni];
    }
    cout << "-------" << endl << endl;
    cout << "median tracking time: " << vTimesTrack[nImages/2] << endl;
    cout << "mean tracking time: " << totaltime/nImages << endl;

    // // Save camera trajectory
    // SLAM.SaveTrajectoryTUM("CameraTrajectory.txt");
    // SLAM.SaveKeyFrameTrajectoryTUM("KeyFrameTrajectory.txt");   

    // return 0;
}

void LoadImages(const string &times, string& file_type, string& basic_rgb, string& basic_depth, vector<string> &vstrImageFilenamesRGB,
                vector<string> &vstrImageFilenamesD, vector<double> &vTimestamps)
{
    string depth_f("/tsukuba_depth_L_");
    string depth_b(".xml");
    string rgb_f("/tsukuba_"); 
    string rgb_m("_L_");
    string rgb_b(".png");

    ifstream time(times.c_str());
    string one_row;
    int sequence = 0;
    double timestamp = 0.0;
    double cycle = 1.0/25.0;
    while(getline(time,one_row)){
        istringstream temp_one_row(one_row);
        string one_ele;
        while(getline(temp_one_row, one_ele, ' ')){
            if(sequence == 0){
                string add_depth = basic_depth + depth_f + one_ele + depth_b;
                string add_rgb = basic_rgb + rgb_f + file_type + rgb_m + one_ele + rgb_b;
                vstrImageFilenamesRGB.push_back(add_rgb);
                vstrImageFilenamesD.push_back(add_depth);
                sequence++;
                // std::cout<<"the address is "<<add_depth.c_str()<<std::endl;
                // cv::Mat imRGB = cv::imread(add_rgb,CV_LOAD_IMAGE_UNCHANGED);

                // cv::namedWindow("show img", cv::WINDOW_AUTOSIZE);
                // cv::imshow("show img", imRGB);

                // cv::waitKey(0);

                // cv::FileStorage fs;
                // fs.open(add_depth, cv::FileStorage::READ);
                // cv::Mat depth_img;
                // fs["depth"]>>depth_img;
                // std::cout<<" the depth is \n"<<depth_img.ptr<float>(0)[0]<<","<<depth_img.ptr<float>(479)[639]<<",row "<<depth_img.rows<<", col "<<depth_img.cols<<std::endl;
                // exit(0);
            }
            else if(sequence == 1){
                vTimestamps.push_back(timestamp);
                timestamp += cycle;
                sequence = 0;
            }
        }
    }
}

