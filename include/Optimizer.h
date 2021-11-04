/**
* This file is part of ORB-SLAM2.
*
* Copyright (C) 2014-2016 Raúl Mur-Artal <raulmur at unizar dot es> (University of Zaragoza)
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

#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "Map.h"
#include "MapPoint.h"
#include "KeyFrame.h"
#include "LoopClosing.h"
#include "Frame.h"
#include "extern_settings.h"

#include "Thirdparty/g2o/g2o/types/types_seven_dof_expmap.h"

namespace ORB_SLAM2
{

class LoopClosing;

class Optimizer
{
public:
    void static BundleAdjustment(const std::vector<KeyFrame*> &vpKF, const std::vector<MapPoint*> &vpMP,
                                 int nIterations = 5, bool *pbStopFlag=NULL, const unsigned long nLoopKF=0,
                                 const bool bRobust = true);
    void static GlobalBundleAdjustemnt(Map* pMap, int nIterations=5, bool *pbStopFlag=NULL,
                                       const unsigned long nLoopKF=0, const bool bRobust = true);
    void static LocalBundleAdjustment(KeyFrame* pKF, bool *pbStopFlag, Map *pMap);
    int static PoseOptimization(Frame* pFrame);

    int static PoseOptimizationNID(Frame *pFrame, Frame *lFrame, const cv::Mat& T_cw);
    int static PoseOptimizationNIDWholeImage(Frame *pFrame, Frame *lFrame, const cv::Mat& T_cw);
    //int static PoseOptimizationIntensity(Frame *pFrame, Frame *lFrame, const cv::Mat& Velocity);
    //float static PoseOptimizationIntensityCoarse(Frame *pFrame, Frame *lFrame, cv::Mat& Velocity);

    bool static PoseOptimizationDirectCoarse(Frame *pFrame, Frame *lFrame, cv::Mat T_cw);
    bool static PoseOptimizationDirect(Frame *pFrame, Frame *lFrame, cv::Mat T_cw);


    // if bFixScale is true, 6DoF optimization (stereo,rgbd), 7DoF otherwise (mono)
    void static OptimizeEssentialGraph(Map* pMap, KeyFrame* pLoopKF, KeyFrame* pCurKF,
                                       const LoopClosing::KeyFrameAndPose &NonCorrectedSim3,
                                       const LoopClosing::KeyFrameAndPose &CorrectedSim3,
                                       const map<KeyFrame *, set<KeyFrame *> > &LoopConnections,
                                       const bool &bFixScale);

    // if bFixScale is true, optimize SE3 (stereo,rgbd), Sim3 otherwise (mono)
    static int OptimizeSim3(KeyFrame* pKF1, KeyFrame* pKF2, std::vector<MapPoint *> &vpMatches1,
                            g2o::Sim3 &g2oS12, const float th2, const bool bFixScale);
private:
    void static Get3dPointAndIntensity(int cell_size, int cell_row_id, int cell_col_id, const cv::Mat& image, const cv::Mat&depth, std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> >& points_3d, Eigen::Matrix<double, Eigen::Dynamic, 1>& intensity, const Eigen::Matrix4d& T_wc, float fx, float fy, float cx, float cy);

    void static Get3dPointAndIntensity(int cell_size, int cell_row_id, int cell_col_id, const cv::Mat& image, const cv::Mat&depth, std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> >& points_3d, Eigen::Matrix<double, Eigen::Dynamic, 1>& intensity, std::vector<cv::KeyPoint>& points_2d,const Eigen::Matrix4d& T_wc, float fx, float fy, float cx, float cy);

};

} //namespace ORB_SLAM

#endif // OPTIMIZER_H
