#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <opencv2/opencv.hpp>
#include <iostream>
#include <stdio.h>

namespace ORB_SLAM2
{
   void ExtractHighGradient(const cv::Mat& im, const cv::Mat& depth, std::vector<cv::KeyPoint>& hg_points, std::vector<float>& hg_depth, const int cell, const int pattern_size);
}