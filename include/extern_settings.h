#pragma once

#include <string.h>
#include <string>
#include <cmath>
#include <vector>

namespace ORB_SLAM2
{
    //b spline coeffcient, b spline derivatives, b spline degree
    extern std::vector<std::vector<std::vector<float> > > bs_;
    extern std::vector<std::vector<std::vector<float> > > bs_der_;
    extern int bs_degree_ ;
    extern int bin_num_ ;
    extern int all_frame_counter_ ;
    extern double all_accu_error_;
    extern int all_available_edge_;
    extern int has_good_guess_counter_;
}// namespace ORB_SLAM2

