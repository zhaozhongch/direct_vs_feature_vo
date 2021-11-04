#include "extern_settings.h"

namespace ORB_SLAM2
{
    //b spline coeffcient, b spline derivatives, b spline degree
    std::vector<std::vector<std::vector<float> > > bs_;
    std::vector<std::vector<std::vector<float> > > bs_der_;
    int bs_degree_ = 3;
    int bin_num_ = 8;
    int all_frame_counter_ = 0;
    double all_accu_error_ = 0.0;
    int all_available_edge_ = 0;
    int has_good_guess_counter_ = 0;
} // namespace ORB_SLAM2
