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

#include "Optimizer.h"

#include "Thirdparty/g2o/g2o/core/block_solver.h"
#include "Thirdparty/g2o/g2o/core/optimization_algorithm_levenberg.h"
#include "Thirdparty/g2o/g2o/solvers/linear_solver_eigen.h"
#include "Thirdparty/g2o/g2o/types/types_six_dof_expmap.h"
#include "Thirdparty/g2o/g2o/core/robust_kernel_impl.h"
#include "Thirdparty/g2o/g2o/solvers/linear_solver_dense.h"
#include "Thirdparty/g2o/g2o/types/types_seven_dof_expmap.h"

#include<Eigen/StdVector>

#include "Converter.h"

#include<mutex>

namespace ORB_SLAM2
{
int one_nid_edge_size_ =500;//should be the same as the SIZE_M in g2o

void Optimizer::GlobalBundleAdjustemnt(Map* pMap, int nIterations, bool* pbStopFlag, const unsigned long nLoopKF, const bool bRobust)
{
    vector<KeyFrame*> vpKFs = pMap->GetAllKeyFrames();
    vector<MapPoint*> vpMP = pMap->GetAllMapPoints();
    BundleAdjustment(vpKFs,vpMP,nIterations,pbStopFlag, nLoopKF, bRobust);
}


void Optimizer::BundleAdjustment(const vector<KeyFrame *> &vpKFs, const vector<MapPoint *> &vpMP,
                                 int nIterations, bool* pbStopFlag, const unsigned long nLoopKF, const bool bRobust)
{
    vector<bool> vbNotIncludedMP;
    vbNotIncludedMP.resize(vpMP.size());

    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    if(pbStopFlag)
        optimizer.setForceStopFlag(pbStopFlag);

    long unsigned int maxKFid = 0;

    // Set KeyFrame vertices
    for(size_t i=0; i<vpKFs.size(); i++)
    {
        KeyFrame* pKF = vpKFs[i];
        if(pKF->isBad())
            continue;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKF->GetPose()));
        vSE3->setId(pKF->mnId);
        vSE3->setFixed(pKF->mnId==0);
        optimizer.addVertex(vSE3);
        if(pKF->mnId>maxKFid)
            maxKFid=pKF->mnId;
    }

    const float thHuber2D = sqrt(5.99);
    const float thHuber3D = sqrt(7.815);

    // Set MapPoint vertices
    for(size_t i=0; i<vpMP.size(); i++)
    {
        MapPoint* pMP = vpMP[i];
        if(pMP->isBad())
            continue;
        g2o::VertexSBAPointXYZ* vPoint = new g2o::VertexSBAPointXYZ();
        vPoint->setEstimate(Converter::toVector3d(pMP->GetWorldPos()));
        const int id = pMP->mnId+maxKFid+1;
        vPoint->setId(id);
        vPoint->setMarginalized(true);
        optimizer.addVertex(vPoint);

       const map<KeyFrame*,size_t> observations = pMP->GetObservations();

        int nEdges = 0;
        //SET EDGES
        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(); mit!=observations.end(); mit++)
        {

            KeyFrame* pKF = mit->first;
            if(pKF->isBad() || pKF->mnId>maxKFid)
                continue;

            nEdges++;

            const cv::KeyPoint &kpUn = pKF->mvKeysUn[mit->second];

            if(pKF->mvuRight[mit->second]<0)
            {
                Eigen::Matrix<double,2,1> obs;
                obs << kpUn.pt.x, kpUn.pt.y;

                g2o::EdgeSE3ProjectXYZ* e = new g2o::EdgeSE3ProjectXYZ();

                e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKF->mnId)));
                e->setMeasurement(obs);
                const float &invSigma2 = pKF->mvInvLevelSigma2[kpUn.octave];
                e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

                if(bRobust)
                {
                    g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuber2D);
                }

                e->fx = pKF->fx;
                e->fy = pKF->fy;
                e->cx = pKF->cx;
                e->cy = pKF->cy;

                optimizer.addEdge(e);
            }
            else
            {
                Eigen::Matrix<double,3,1> obs;
                const float kp_ur = pKF->mvuRight[mit->second];
                obs << kpUn.pt.x, kpUn.pt.y, kp_ur;

                g2o::EdgeStereoSE3ProjectXYZ* e = new g2o::EdgeStereoSE3ProjectXYZ();

                e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKF->mnId)));
                e->setMeasurement(obs);
                const float &invSigma2 = pKF->mvInvLevelSigma2[kpUn.octave];
                Eigen::Matrix3d Info = Eigen::Matrix3d::Identity()*invSigma2;
                e->setInformation(Info);

                if(bRobust)
                {
                    g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuber3D);
                }

                e->fx = pKF->fx;
                e->fy = pKF->fy;
                e->cx = pKF->cx;
                e->cy = pKF->cy;
                e->bf = pKF->mbf;

                optimizer.addEdge(e);
            }
        }

        if(nEdges==0)
        {
            optimizer.removeVertex(vPoint);
            vbNotIncludedMP[i]=true;
        }
        else
        {
            vbNotIncludedMP[i]=false;
        }
    }

    // Optimize!
    optimizer.initializeOptimization();
    optimizer.optimize(nIterations);

    // Recover optimized data

    //Keyframes
    for(size_t i=0; i<vpKFs.size(); i++)
    {
        KeyFrame* pKF = vpKFs[i];
        if(pKF->isBad())
            continue;
        g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(pKF->mnId));
        g2o::SE3Quat SE3quat = vSE3->estimate();
        if(nLoopKF==0)
        {
            pKF->SetPose(Converter::toCvMat(SE3quat));
        }
        else
        {
            pKF->mTcwGBA.create(4,4,CV_32F);
            Converter::toCvMat(SE3quat).copyTo(pKF->mTcwGBA);
            pKF->mnBAGlobalForKF = nLoopKF;
        }
    }

    //Points
    for(size_t i=0; i<vpMP.size(); i++)
    {
        if(vbNotIncludedMP[i])
            continue;

        MapPoint* pMP = vpMP[i];

        if(pMP->isBad())
            continue;
        g2o::VertexSBAPointXYZ* vPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(pMP->mnId+maxKFid+1));

        if(nLoopKF==0)
        {
            pMP->SetWorldPos(Converter::toCvMat(vPoint->estimate()));
            pMP->UpdateNormalAndDepth();
        }
        else
        {
            pMP->mPosGBA.create(3,1,CV_32F);
            Converter::toCvMat(vPoint->estimate()).copyTo(pMP->mPosGBA);
            pMP->mnBAGlobalForKF = nLoopKF;
        }
    }

}

int Optimizer::PoseOptimization(Frame *pFrame)
{
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    //optimizer.setVerbose(true);

    int nInitialCorrespondences=0;

    // Set Frame vertex
    g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
    //std::cout<<"before orbslam optimization the matrix in world frame is \n"<<pFrame->mTcw<<std::endl;
    vSE3->setEstimate(Converter::toSE3Quat(pFrame->mTcw));
    vSE3->setId(0);
    vSE3->setFixed(false);

    optimizer.addVertex(vSE3);

    // Set MapPoint vertices
    const int N = pFrame->N;

    vector<g2o::EdgeSE3ProjectXYZOnlyPose*> vpEdgesMono;
    vector<size_t> vnIndexEdgeMono;
    vpEdgesMono.reserve(N);
    vnIndexEdgeMono.reserve(N);

    vector<g2o::EdgeStereoSE3ProjectXYZOnlyPose*> vpEdgesStereo;
    vector<size_t> vnIndexEdgeStereo;
    vpEdgesStereo.reserve(N);
    vnIndexEdgeStereo.reserve(N);

    const float deltaMono = sqrt(5.991);
    const float deltaStereo = sqrt(7.815);


    {
    unique_lock<mutex> lock(MapPoint::mGlobalMutex);

    for(int i=0; i<N; i++)
    {
        MapPoint* pMP = pFrame->mvpMapPoints[i];
        if(pMP)
        {
            // Monocular observation
            if(pFrame->mvuRight[i]<0)
            {
                nInitialCorrespondences++;
                pFrame->mvbOutlier[i] = false;

                Eigen::Matrix<double,2,1> obs;
                const cv::KeyPoint &kpUn = pFrame->mvKeysUn[i];
                obs << kpUn.pt.x, kpUn.pt.y;

                g2o::EdgeSE3ProjectXYZOnlyPose* e = new g2o::EdgeSE3ProjectXYZOnlyPose();

                e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
                e->setMeasurement(obs);
                const float invSigma2 = pFrame->mvInvLevelSigma2[kpUn.octave];
                e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

                g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                e->setRobustKernel(rk);
                rk->setDelta(deltaMono);

                e->fx = pFrame->fx;
                e->fy = pFrame->fy;
                e->cx = pFrame->cx;
                e->cy = pFrame->cy;
                cv::Mat Xw = pMP->GetWorldPos();
                e->Xw[0] = Xw.at<float>(0);
                e->Xw[1] = Xw.at<float>(1);
                e->Xw[2] = Xw.at<float>(2);

                optimizer.addEdge(e);

                vpEdgesMono.push_back(e);
                vnIndexEdgeMono.push_back(i);
            }
            else  // Stereo observation
            {
                nInitialCorrespondences++;
                pFrame->mvbOutlier[i] = false;

                //SET EDGE
                Eigen::Matrix<double,3,1> obs;
                const cv::KeyPoint &kpUn = pFrame->mvKeysUn[i];
                const float &kp_ur = pFrame->mvuRight[i];
                obs << kpUn.pt.x, kpUn.pt.y, kp_ur;

                g2o::EdgeStereoSE3ProjectXYZOnlyPose* e = new g2o::EdgeStereoSE3ProjectXYZOnlyPose();

                e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
                e->setMeasurement(obs);
                const float invSigma2 = pFrame->mvInvLevelSigma2[kpUn.octave];
                Eigen::Matrix3d Info = Eigen::Matrix3d::Identity()*invSigma2;
                e->setInformation(Info);

                g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                e->setRobustKernel(rk);
                rk->setDelta(deltaStereo);

                e->fx = pFrame->fx;
                e->fy = pFrame->fy;
                e->cx = pFrame->cx;
                e->cy = pFrame->cy;
                e->bf = pFrame->mbf;
                cv::Mat Xw = pMP->GetWorldPos();
                e->Xw[0] = Xw.at<float>(0);
                e->Xw[1] = Xw.at<float>(1);
                e->Xw[2] = Xw.at<float>(2);

                optimizer.addEdge(e);

                vpEdgesStereo.push_back(e);
                vnIndexEdgeStereo.push_back(i);
            }
        }

    }
    }


    if(nInitialCorrespondences<3)
        return 0;

    // We perform 4 optimizations, after each optimization we classify observation as inlier/outlier
    // At the next optimization, outliers are not included, but at the end they can be classified as inliers again.
    const float chi2Mono[4]={5.991,5.991,5.991,5.991};
    const float chi2Stereo[4]={7.815,7.815,7.815, 7.815};
    const int its[4]={10,10,10,10};

    int nBad=0;
    for(size_t it=0; it<4; it++)
    {

        vSE3->setEstimate(Converter::toSE3Quat(pFrame->mTcw));
        optimizer.initializeOptimization(0);
        //std::cout<<"enter optimization...................................... "<< it<<std::endl;
        optimizer.optimize(its[it]);
        //std::cout<<"finish optimization...................................... "<< it<<std::endl;

        nBad=0;
        for(size_t i=0, iend=vpEdgesMono.size(); i<iend; i++)
        {
            g2o::EdgeSE3ProjectXYZOnlyPose* e = vpEdgesMono[i];

            const size_t idx = vnIndexEdgeMono[i];

            if(pFrame->mvbOutlier[idx])
            {
                e->computeError();
            }

            const float chi2 = e->chi2();

            if(chi2>chi2Mono[it])
            {                
                pFrame->mvbOutlier[idx]=true;
                e->setLevel(1);
                nBad++;
            }
            else
            {
                pFrame->mvbOutlier[idx]=false;
                e->setLevel(0);
            }

            if(it==2)
                e->setRobustKernel(0);
        }

        for(size_t i=0, iend=vpEdgesStereo.size(); i<iend; i++)
        {
            g2o::EdgeStereoSE3ProjectXYZOnlyPose* e = vpEdgesStereo[i];

            const size_t idx = vnIndexEdgeStereo[i];

            if(pFrame->mvbOutlier[idx])
            {
                e->computeError();
            }

            const float chi2 = e->chi2();

            if(chi2>chi2Stereo[it])
            {
                pFrame->mvbOutlier[idx]=true;
                e->setLevel(1);
                nBad++;
            }
            else
            {                
                e->setLevel(0);
                pFrame->mvbOutlier[idx]=false;
            }

            if(it==2)
                e->setRobustKernel(0);
        }

        if(optimizer.edges().size()<10)
            break;
    }    

    //std::cout<<"\n \n"<<std::endl;

    // Recover optimized pose and return number of inliers
    g2o::VertexSE3Expmap* vSE3_recov = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(0));
    g2o::SE3Quat SE3quat_recov = vSE3_recov->estimate();
    cv::Mat pose = Converter::toCvMat(SE3quat_recov);
    pFrame->SetPose(pose);

    // std::cout<<"after orbslam optimization the matrix in world frame is \n"<<pose<<std::endl;
    // std::cout<<std::endl;

    return nInitialCorrespondences-nBad;
}

void Optimizer::LocalBundleAdjustment(KeyFrame *pKF, bool* pbStopFlag, Map* pMap)
{    
    // Local KeyFrames: First Breath Search from Current Keyframe
    list<KeyFrame*> lLocalKeyFrames;

    lLocalKeyFrames.push_back(pKF);
    pKF->mnBALocalForKF = pKF->mnId;

    const vector<KeyFrame*> vNeighKFs = pKF->GetVectorCovisibleKeyFrames();
    for(int i=0, iend=vNeighKFs.size(); i<iend; i++)
    {
        KeyFrame* pKFi = vNeighKFs[i];
        pKFi->mnBALocalForKF = pKF->mnId;
        if(!pKFi->isBad())
            lLocalKeyFrames.push_back(pKFi);
    }

    // Local MapPoints seen in Local KeyFrames
    list<MapPoint*> lLocalMapPoints;
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin() , lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        vector<MapPoint*> vpMPs = (*lit)->GetMapPointMatches();
        for(vector<MapPoint*>::iterator vit=vpMPs.begin(), vend=vpMPs.end(); vit!=vend; vit++)
        {
            MapPoint* pMP = *vit;
            if(pMP)
                if(!pMP->isBad())
                    if(pMP->mnBALocalForKF!=pKF->mnId)
                    {
                        lLocalMapPoints.push_back(pMP);
                        pMP->mnBALocalForKF=pKF->mnId;
                    }
        }
    }

    // Fixed Keyframes. Keyframes that see Local MapPoints but that are not Local Keyframes
    list<KeyFrame*> lFixedCameras;
    for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        map<KeyFrame*,size_t> observations = (*lit)->GetObservations();
        for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(pKFi->mnBALocalForKF!=pKF->mnId && pKFi->mnBAFixedForKF!=pKF->mnId)
            {                
                pKFi->mnBAFixedForKF=pKF->mnId;
                if(!pKFi->isBad())
                    lFixedCameras.push_back(pKFi);
            }
        }
    }

    // Setup optimizer
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    //optimizer.setVerbose(true);
    double opt_time = (double)cv::getTickCount();

    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    if(pbStopFlag)
        optimizer.setForceStopFlag(pbStopFlag);

    unsigned long maxKFid = 0;

    //count vertex and edges
    //int vcount = 0;
    //int ecount = 0;

    // Set Local KeyFrame vertices
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        //vcount+=6;
        KeyFrame* pKFi = *lit;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
        vSE3->setId(pKFi->mnId);
        vSE3->setFixed(pKFi->mnId==0);
        optimizer.addVertex(vSE3);
        if(pKFi->mnId>maxKFid)
            maxKFid=pKFi->mnId;
    }

    // Set Fixed KeyFrame vertices
    for(list<KeyFrame*>::iterator lit=lFixedCameras.begin(), lend=lFixedCameras.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
        vSE3->setId(pKFi->mnId);
        vSE3->setFixed(true);
        optimizer.addVertex(vSE3);
        if(pKFi->mnId>maxKFid)
            maxKFid=pKFi->mnId;
    }

    // Set MapPoint vertices
    const int nExpectedSize = (lLocalKeyFrames.size()+lFixedCameras.size())*lLocalMapPoints.size();

    vector<g2o::EdgeSE3ProjectXYZ*> vpEdgesMono;
    vpEdgesMono.reserve(nExpectedSize);

    vector<KeyFrame*> vpEdgeKFMono;
    vpEdgeKFMono.reserve(nExpectedSize);

    vector<MapPoint*> vpMapPointEdgeMono;
    vpMapPointEdgeMono.reserve(nExpectedSize);

    vector<g2o::EdgeStereoSE3ProjectXYZ*> vpEdgesStereo;
    vpEdgesStereo.reserve(nExpectedSize);

    vector<KeyFrame*> vpEdgeKFStereo;
    vpEdgeKFStereo.reserve(nExpectedSize);

    vector<MapPoint*> vpMapPointEdgeStereo;
    vpMapPointEdgeStereo.reserve(nExpectedSize);

    const float thHuberMono = sqrt(5.991);
    const float thHuberStereo = sqrt(7.815);

    for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        MapPoint* pMP = *lit;
        g2o::VertexSBAPointXYZ* vPoint = new g2o::VertexSBAPointXYZ();
        vPoint->setEstimate(Converter::toVector3d(pMP->GetWorldPos()));
        int id = pMP->mnId+maxKFid+1;
        vPoint->setId(id);
        vPoint->setMarginalized(true);
        optimizer.addVertex(vPoint);
        //vcount+=3;

        const map<KeyFrame*,size_t> observations = pMP->GetObservations();

        //Set edges
        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(!pKFi->isBad())
            {                
                const cv::KeyPoint &kpUn = pKFi->mvKeysUn[mit->second];

                // Monocular observation
                if(pKFi->mvuRight[mit->second]<0)
                {
                    Eigen::Matrix<double,2,1> obs;
                    obs << kpUn.pt.x, kpUn.pt.y;

                    g2o::EdgeSE3ProjectXYZ* e = new g2o::EdgeSE3ProjectXYZ();

                    e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                    e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
                    e->setMeasurement(obs);
                    const float &invSigma2 = pKFi->mvInvLevelSigma2[kpUn.octave];
                    e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

                    g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuberMono);

                    e->fx = pKFi->fx;
                    e->fy = pKFi->fy;
                    e->cx = pKFi->cx;
                    e->cy = pKFi->cy;

                    optimizer.addEdge(e);
                    vpEdgesMono.push_back(e);
                    vpEdgeKFMono.push_back(pKFi);
                    vpMapPointEdgeMono.push_back(pMP);
                }
                else // Stereo observation
                {
                    Eigen::Matrix<double,3,1> obs;
                    const float kp_ur = pKFi->mvuRight[mit->second];
                    obs << kpUn.pt.x, kpUn.pt.y, kp_ur;

                    g2o::EdgeStereoSE3ProjectXYZ* e = new g2o::EdgeStereoSE3ProjectXYZ();

                    e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                    e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
                    e->setMeasurement(obs);
                    const float &invSigma2 = pKFi->mvInvLevelSigma2[kpUn.octave];
                    Eigen::Matrix3d Info = Eigen::Matrix3d::Identity()*invSigma2;
                    e->setInformation(Info);

                    g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuberStereo);

                    e->fx = pKFi->fx;
                    e->fy = pKFi->fy;
                    e->cx = pKFi->cx;
                    e->cy = pKFi->cy;
                    e->bf = pKFi->mbf;

                    optimizer.addEdge(e);
                    vpEdgesStereo.push_back(e);
                    vpEdgeKFStereo.push_back(pKFi);
                    vpMapPointEdgeStereo.push_back(pMP);
                }
            }
        }
    }


    if(pbStopFlag)
        if(*pbStopFlag)
            return;

    optimizer.initializeOptimization();
    optimizer.optimize(5);

    bool bDoMore= true;

    if(pbStopFlag)
        if(*pbStopFlag)
            bDoMore = false;

    if(bDoMore)
    {

    // Check inlier observations
    for(size_t i=0, iend=vpEdgesMono.size(); i<iend;i++)
    {
        g2o::EdgeSE3ProjectXYZ* e = vpEdgesMono[i];
        MapPoint* pMP = vpMapPointEdgeMono[i];

        if(pMP->isBad())
            continue;

        if(e->chi2()>5.991 || !e->isDepthPositive())
        {
            e->setLevel(1);
        }

        e->setRobustKernel(0);
    }

    for(size_t i=0, iend=vpEdgesStereo.size(); i<iend;i++)
    {
        g2o::EdgeStereoSE3ProjectXYZ* e = vpEdgesStereo[i];
        MapPoint* pMP = vpMapPointEdgeStereo[i];

        if(pMP->isBad())
            continue;

        if(e->chi2()>7.815 || !e->isDepthPositive())
        {
            e->setLevel(1);
        }

        e->setRobustKernel(0);
    }

    // Optimize again without the outliers

    optimizer.initializeOptimization(0);
    optimizer.optimize(10);

    }

    opt_time = ((double)cv::getTickCount() - opt_time)/cv::getTickFrequency();

    //std::cout<<"in localBA we have ........"<<vcount<<".....parameters to optimize, pose number is  "<<lLocalKeyFrames.size()<<", point number is "<<(vcount - lLocalKeyFrames.size()*6)/3<<std::endl;
    //std::cout<<"local BA optimization time is................. "<<opt_time<<std::endl;

    vector<pair<KeyFrame*,MapPoint*> > vToErase;
    vToErase.reserve(vpEdgesMono.size()+vpEdgesStereo.size());

    // Check inlier observations       
    for(size_t i=0, iend=vpEdgesMono.size(); i<iend;i++)
    {
        g2o::EdgeSE3ProjectXYZ* e = vpEdgesMono[i];
        MapPoint* pMP = vpMapPointEdgeMono[i];

        if(pMP->isBad())
            continue;

        if(e->chi2()>5.991 || !e->isDepthPositive())
        {
            KeyFrame* pKFi = vpEdgeKFMono[i];
            vToErase.push_back(make_pair(pKFi,pMP));
        }
    }

    for(size_t i=0, iend=vpEdgesStereo.size(); i<iend;i++)
    {
        g2o::EdgeStereoSE3ProjectXYZ* e = vpEdgesStereo[i];
        MapPoint* pMP = vpMapPointEdgeStereo[i];

        if(pMP->isBad())
            continue;

        if(e->chi2()>7.815 || !e->isDepthPositive())
        {
            KeyFrame* pKFi = vpEdgeKFStereo[i];
            vToErase.push_back(make_pair(pKFi,pMP));
        }
    }

    // Get Map Mutex
    unique_lock<mutex> lock(pMap->mMutexMapUpdate);

    if(!vToErase.empty())
    {
        for(size_t i=0;i<vToErase.size();i++)
        {
            KeyFrame* pKFi = vToErase[i].first;
            MapPoint* pMPi = vToErase[i].second;
            pKFi->EraseMapPointMatch(pMPi);
            pMPi->EraseObservation(pKFi);
        }
    }

    // Recover optimized data

    //Keyframes
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        KeyFrame* pKF = *lit;
        g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(pKF->mnId));
        g2o::SE3Quat SE3quat = vSE3->estimate();
        pKF->SetPose(Converter::toCvMat(SE3quat));
    }

    //Points
    for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        MapPoint* pMP = *lit;
        g2o::VertexSBAPointXYZ* vPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(pMP->mnId+maxKFid+1));
        pMP->SetWorldPos(Converter::toCvMat(vPoint->estimate()));
        pMP->UpdateNormalAndDepth();
    }
}


void Optimizer::OptimizeEssentialGraph(Map* pMap, KeyFrame* pLoopKF, KeyFrame* pCurKF,
                                       const LoopClosing::KeyFrameAndPose &NonCorrectedSim3,
                                       const LoopClosing::KeyFrameAndPose &CorrectedSim3,
                                       const map<KeyFrame *, set<KeyFrame *> > &LoopConnections, const bool &bFixScale)
{
    // Setup optimizer
    g2o::SparseOptimizer optimizer;
    optimizer.setVerbose(false);
    g2o::BlockSolver_7_3::LinearSolverType * linearSolver =
           new g2o::LinearSolverEigen<g2o::BlockSolver_7_3::PoseMatrixType>();
    g2o::BlockSolver_7_3 * solver_ptr= new g2o::BlockSolver_7_3(linearSolver);
    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);

    solver->setUserLambdaInit(1e-16);
    optimizer.setAlgorithm(solver);

    const vector<KeyFrame*> vpKFs = pMap->GetAllKeyFrames();
    const vector<MapPoint*> vpMPs = pMap->GetAllMapPoints();

    const unsigned int nMaxKFid = pMap->GetMaxKFid();

    vector<g2o::Sim3,Eigen::aligned_allocator<g2o::Sim3> > vScw(nMaxKFid+1);
    vector<g2o::Sim3,Eigen::aligned_allocator<g2o::Sim3> > vCorrectedSwc(nMaxKFid+1);
    vector<g2o::VertexSim3Expmap*> vpVertices(nMaxKFid+1);

    const int minFeat = 100;

    // Set KeyFrame vertices
    for(size_t i=0, iend=vpKFs.size(); i<iend;i++)
    {
        KeyFrame* pKF = vpKFs[i];
        if(pKF->isBad())
            continue;
        g2o::VertexSim3Expmap* VSim3 = new g2o::VertexSim3Expmap();

        const int nIDi = pKF->mnId;

        LoopClosing::KeyFrameAndPose::const_iterator it = CorrectedSim3.find(pKF);

        if(it!=CorrectedSim3.end())
        {
            vScw[nIDi] = it->second;
            VSim3->setEstimate(it->second);
        }
        else
        {
            Eigen::Matrix<double,3,3> Rcw = Converter::toMatrix3d(pKF->GetRotation());
            Eigen::Matrix<double,3,1> tcw = Converter::toVector3d(pKF->GetTranslation());
            g2o::Sim3 Siw(Rcw,tcw,1.0);
            vScw[nIDi] = Siw;
            VSim3->setEstimate(Siw);
        }

        if(pKF==pLoopKF)
            VSim3->setFixed(true);

        VSim3->setId(nIDi);
        VSim3->setMarginalized(false);
        VSim3->_fix_scale = bFixScale;

        optimizer.addVertex(VSim3);

        vpVertices[nIDi]=VSim3;
    }


    set<pair<long unsigned int,long unsigned int> > sInsertedEdges;

    const Eigen::Matrix<double,7,7> matLambda = Eigen::Matrix<double,7,7>::Identity();

    // Set Loop edges
    for(map<KeyFrame *, set<KeyFrame *> >::const_iterator mit = LoopConnections.begin(), mend=LoopConnections.end(); mit!=mend; mit++)
    {
        KeyFrame* pKF = mit->first;
        const long unsigned int nIDi = pKF->mnId;
        const set<KeyFrame*> &spConnections = mit->second;
        const g2o::Sim3 Siw = vScw[nIDi];
        const g2o::Sim3 Swi = Siw.inverse();

        for(set<KeyFrame*>::const_iterator sit=spConnections.begin(), send=spConnections.end(); sit!=send; sit++)
        {
            const long unsigned int nIDj = (*sit)->mnId;
            if((nIDi!=pCurKF->mnId || nIDj!=pLoopKF->mnId) && pKF->GetWeight(*sit)<minFeat)
                continue;

            const g2o::Sim3 Sjw = vScw[nIDj];
            const g2o::Sim3 Sji = Sjw * Swi;

            g2o::EdgeSim3* e = new g2o::EdgeSim3();
            e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDj)));
            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDi)));
            e->setMeasurement(Sji);

            e->information() = matLambda;

            optimizer.addEdge(e);

            sInsertedEdges.insert(make_pair(min(nIDi,nIDj),max(nIDi,nIDj)));
        }
    }

    // Set normal edges
    for(size_t i=0, iend=vpKFs.size(); i<iend; i++)
    {
        KeyFrame* pKF = vpKFs[i];

        const int nIDi = pKF->mnId;

        g2o::Sim3 Swi;

        LoopClosing::KeyFrameAndPose::const_iterator iti = NonCorrectedSim3.find(pKF);

        if(iti!=NonCorrectedSim3.end())
            Swi = (iti->second).inverse();
        else
            Swi = vScw[nIDi].inverse();

        KeyFrame* pParentKF = pKF->GetParent();

        // Spanning tree edge
        if(pParentKF)
        {
            int nIDj = pParentKF->mnId;

            g2o::Sim3 Sjw;

            LoopClosing::KeyFrameAndPose::const_iterator itj = NonCorrectedSim3.find(pParentKF);

            if(itj!=NonCorrectedSim3.end())
                Sjw = itj->second;
            else
                Sjw = vScw[nIDj];

            g2o::Sim3 Sji = Sjw * Swi;

            g2o::EdgeSim3* e = new g2o::EdgeSim3();
            e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDj)));
            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDi)));
            e->setMeasurement(Sji);

            e->information() = matLambda;
            optimizer.addEdge(e);
        }

        // Loop edges
        const set<KeyFrame*> sLoopEdges = pKF->GetLoopEdges();
        for(set<KeyFrame*>::const_iterator sit=sLoopEdges.begin(), send=sLoopEdges.end(); sit!=send; sit++)
        {
            KeyFrame* pLKF = *sit;
            if(pLKF->mnId<pKF->mnId)
            {
                g2o::Sim3 Slw;

                LoopClosing::KeyFrameAndPose::const_iterator itl = NonCorrectedSim3.find(pLKF);

                if(itl!=NonCorrectedSim3.end())
                    Slw = itl->second;
                else
                    Slw = vScw[pLKF->mnId];

                g2o::Sim3 Sli = Slw * Swi;
                g2o::EdgeSim3* el = new g2o::EdgeSim3();
                el->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pLKF->mnId)));
                el->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDi)));
                el->setMeasurement(Sli);
                el->information() = matLambda;
                optimizer.addEdge(el);
            }
        }

        // Covisibility graph edges
        const vector<KeyFrame*> vpConnectedKFs = pKF->GetCovisiblesByWeight(minFeat);
        for(vector<KeyFrame*>::const_iterator vit=vpConnectedKFs.begin(); vit!=vpConnectedKFs.end(); vit++)
        {
            KeyFrame* pKFn = *vit;
            if(pKFn && pKFn!=pParentKF && !pKF->hasChild(pKFn) && !sLoopEdges.count(pKFn))
            {
                if(!pKFn->isBad() && pKFn->mnId<pKF->mnId)
                {
                    if(sInsertedEdges.count(make_pair(min(pKF->mnId,pKFn->mnId),max(pKF->mnId,pKFn->mnId))))
                        continue;

                    g2o::Sim3 Snw;

                    LoopClosing::KeyFrameAndPose::const_iterator itn = NonCorrectedSim3.find(pKFn);

                    if(itn!=NonCorrectedSim3.end())
                        Snw = itn->second;
                    else
                        Snw = vScw[pKFn->mnId];

                    g2o::Sim3 Sni = Snw * Swi;

                    g2o::EdgeSim3* en = new g2o::EdgeSim3();
                    en->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFn->mnId)));
                    en->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDi)));
                    en->setMeasurement(Sni);
                    en->information() = matLambda;
                    optimizer.addEdge(en);
                }
            }
        }
    }

    // Optimize!
    optimizer.initializeOptimization();
    optimizer.optimize(20);

    unique_lock<mutex> lock(pMap->mMutexMapUpdate);

    // SE3 Pose Recovering. Sim3:[sR t;0 1] -> SE3:[R t/s;0 1]
    for(size_t i=0;i<vpKFs.size();i++)
    {
        KeyFrame* pKFi = vpKFs[i];

        const int nIDi = pKFi->mnId;

        g2o::VertexSim3Expmap* VSim3 = static_cast<g2o::VertexSim3Expmap*>(optimizer.vertex(nIDi));
        g2o::Sim3 CorrectedSiw =  VSim3->estimate();
        vCorrectedSwc[nIDi]=CorrectedSiw.inverse();
        Eigen::Matrix3d eigR = CorrectedSiw.rotation().toRotationMatrix();
        Eigen::Vector3d eigt = CorrectedSiw.translation();
        double s = CorrectedSiw.scale();

        eigt *=(1./s); //[R t/s;0 1]

        cv::Mat Tiw = Converter::toCvSE3(eigR,eigt);

        pKFi->SetPose(Tiw);
    }

    // Correct points. Transform to "non-optimized" reference keyframe pose and transform back with optimized pose
    for(size_t i=0, iend=vpMPs.size(); i<iend; i++)
    {
        MapPoint* pMP = vpMPs[i];

        if(pMP->isBad())
            continue;

        int nIDr;
        if(pMP->mnCorrectedByKF==pCurKF->mnId)
        {
            nIDr = pMP->mnCorrectedReference;
        }
        else
        {
            KeyFrame* pRefKF = pMP->GetReferenceKeyFrame();
            nIDr = pRefKF->mnId;
        }


        g2o::Sim3 Srw = vScw[nIDr];
        g2o::Sim3 correctedSwr = vCorrectedSwc[nIDr];

        cv::Mat P3Dw = pMP->GetWorldPos();
        Eigen::Matrix<double,3,1> eigP3Dw = Converter::toVector3d(P3Dw);
        Eigen::Matrix<double,3,1> eigCorrectedP3Dw = correctedSwr.map(Srw.map(eigP3Dw));

        cv::Mat cvCorrectedP3Dw = Converter::toCvMat(eigCorrectedP3Dw);
        pMP->SetWorldPos(cvCorrectedP3Dw);

        pMP->UpdateNormalAndDepth();
    }
}

int Optimizer::OptimizeSim3(KeyFrame *pKF1, KeyFrame *pKF2, vector<MapPoint *> &vpMatches1, g2o::Sim3 &g2oS12, const float th2, const bool bFixScale)
{
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolverX::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolverX::PoseMatrixType>();

    g2o::BlockSolverX * solver_ptr = new g2o::BlockSolverX(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    // Calibration
    const cv::Mat &K1 = pKF1->mK;
    const cv::Mat &K2 = pKF2->mK;

    // Camera poses
    const cv::Mat R1w = pKF1->GetRotation();
    const cv::Mat t1w = pKF1->GetTranslation();
    const cv::Mat R2w = pKF2->GetRotation();
    const cv::Mat t2w = pKF2->GetTranslation();

    // Set Sim3 vertex
    g2o::VertexSim3Expmap * vSim3 = new g2o::VertexSim3Expmap();    
    vSim3->_fix_scale=bFixScale;
    vSim3->setEstimate(g2oS12);
    vSim3->setId(0);
    vSim3->setFixed(false);
    vSim3->_principle_point1[0] = K1.at<float>(0,2);
    vSim3->_principle_point1[1] = K1.at<float>(1,2);
    vSim3->_focal_length1[0] = K1.at<float>(0,0);
    vSim3->_focal_length1[1] = K1.at<float>(1,1);
    vSim3->_principle_point2[0] = K2.at<float>(0,2);
    vSim3->_principle_point2[1] = K2.at<float>(1,2);
    vSim3->_focal_length2[0] = K2.at<float>(0,0);
    vSim3->_focal_length2[1] = K2.at<float>(1,1);
    optimizer.addVertex(vSim3);

    // Set MapPoint vertices
    const int N = vpMatches1.size();
    const vector<MapPoint*> vpMapPoints1 = pKF1->GetMapPointMatches();
    vector<g2o::EdgeSim3ProjectXYZ*> vpEdges12;
    vector<g2o::EdgeInverseSim3ProjectXYZ*> vpEdges21;
    vector<size_t> vnIndexEdge;

    vnIndexEdge.reserve(2*N);
    vpEdges12.reserve(2*N);
    vpEdges21.reserve(2*N);

    const float deltaHuber = sqrt(th2);

    int nCorrespondences = 0;

    for(int i=0; i<N; i++)
    {
        if(!vpMatches1[i])
            continue;

        MapPoint* pMP1 = vpMapPoints1[i];
        MapPoint* pMP2 = vpMatches1[i];

        const int id1 = 2*i+1;
        const int id2 = 2*(i+1);

        const int i2 = pMP2->GetIndexInKeyFrame(pKF2);

        if(pMP1 && pMP2)
        {
            if(!pMP1->isBad() && !pMP2->isBad() && i2>=0)
            {
                g2o::VertexSBAPointXYZ* vPoint1 = new g2o::VertexSBAPointXYZ();
                cv::Mat P3D1w = pMP1->GetWorldPos();
                cv::Mat P3D1c = R1w*P3D1w + t1w;
                vPoint1->setEstimate(Converter::toVector3d(P3D1c));
                vPoint1->setId(id1);
                vPoint1->setFixed(true);
                optimizer.addVertex(vPoint1);

                g2o::VertexSBAPointXYZ* vPoint2 = new g2o::VertexSBAPointXYZ();
                cv::Mat P3D2w = pMP2->GetWorldPos();
                cv::Mat P3D2c = R2w*P3D2w + t2w;
                vPoint2->setEstimate(Converter::toVector3d(P3D2c));
                vPoint2->setId(id2);
                vPoint2->setFixed(true);
                optimizer.addVertex(vPoint2);
            }
            else
                continue;
        }
        else
            continue;

        nCorrespondences++;

        // Set edge x1 = S12*X2
        Eigen::Matrix<double,2,1> obs1;
        const cv::KeyPoint &kpUn1 = pKF1->mvKeysUn[i];
        obs1 << kpUn1.pt.x, kpUn1.pt.y;

        g2o::EdgeSim3ProjectXYZ* e12 = new g2o::EdgeSim3ProjectXYZ();
        e12->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id2)));
        e12->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
        e12->setMeasurement(obs1);
        const float &invSigmaSquare1 = pKF1->mvInvLevelSigma2[kpUn1.octave];
        e12->setInformation(Eigen::Matrix2d::Identity()*invSigmaSquare1);

        g2o::RobustKernelHuber* rk1 = new g2o::RobustKernelHuber;
        e12->setRobustKernel(rk1);
        rk1->setDelta(deltaHuber);
        optimizer.addEdge(e12);

        // Set edge x2 = S21*X1
        Eigen::Matrix<double,2,1> obs2;
        const cv::KeyPoint &kpUn2 = pKF2->mvKeysUn[i2];
        obs2 << kpUn2.pt.x, kpUn2.pt.y;

        g2o::EdgeInverseSim3ProjectXYZ* e21 = new g2o::EdgeInverseSim3ProjectXYZ();

        e21->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id1)));
        e21->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
        e21->setMeasurement(obs2);
        float invSigmaSquare2 = pKF2->mvInvLevelSigma2[kpUn2.octave];
        e21->setInformation(Eigen::Matrix2d::Identity()*invSigmaSquare2);

        g2o::RobustKernelHuber* rk2 = new g2o::RobustKernelHuber;
        e21->setRobustKernel(rk2);
        rk2->setDelta(deltaHuber);
        optimizer.addEdge(e21);

        vpEdges12.push_back(e12);
        vpEdges21.push_back(e21);
        vnIndexEdge.push_back(i);
    }

    // Optimize!
    optimizer.initializeOptimization();
    optimizer.optimize(5);

    // Check inliers
    int nBad=0;
    for(size_t i=0; i<vpEdges12.size();i++)
    {
        g2o::EdgeSim3ProjectXYZ* e12 = vpEdges12[i];
        g2o::EdgeInverseSim3ProjectXYZ* e21 = vpEdges21[i];
        if(!e12 || !e21)
            continue;

        if(e12->chi2()>th2 || e21->chi2()>th2)
        {
            size_t idx = vnIndexEdge[i];
            vpMatches1[idx]=static_cast<MapPoint*>(NULL);
            optimizer.removeEdge(e12);
            optimizer.removeEdge(e21);
            vpEdges12[i]=static_cast<g2o::EdgeSim3ProjectXYZ*>(NULL);
            vpEdges21[i]=static_cast<g2o::EdgeInverseSim3ProjectXYZ*>(NULL);
            nBad++;
        }
    }

    int nMoreIterations;
    if(nBad>0)
        nMoreIterations=10;
    else
        nMoreIterations=5;

    if(nCorrespondences-nBad<10)
        return 0;

    // Optimize again only with inliers

    optimizer.initializeOptimization();
    optimizer.optimize(nMoreIterations);

    int nIn = 0;
    for(size_t i=0; i<vpEdges12.size();i++)
    {
        g2o::EdgeSim3ProjectXYZ* e12 = vpEdges12[i];
        g2o::EdgeInverseSim3ProjectXYZ* e21 = vpEdges21[i];
        if(!e12 || !e21)
            continue;

        if(e12->chi2()>th2 || e21->chi2()>th2)
        {
            size_t idx = vnIndexEdge[i];
            vpMatches1[idx]=static_cast<MapPoint*>(NULL);
        }
        else
            nIn++;
    }

    // Recover optimized Sim3
    g2o::VertexSim3Expmap* vSim3_recov = static_cast<g2o::VertexSim3Expmap*>(optimizer.vertex(0));
    g2oS12= vSim3_recov->estimate();

    return nIn;
}

//measurement should be all 2d points in one cell. Compute the NID in the g2o function
int Optimizer::PoseOptimizationNID(Frame *pFrame, Frame *lFrame, const cv::Mat& T_cw)
{
    int level = 1;
    int cell_size = pFrame->cell_  * pFrame->cell_ / ((level + 1) * (level + 1));
    
    Eigen::MatrixXd hg_intensity_last(one_nid_edge_size_, cell_size);
    //std::vector<std::vector<float> > hg_intensity_last(cell_size);
    std::vector<std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > > hg_keys_3d_last(cell_size);
    std::vector<int> cell_counter(cell_size, 0);

    Eigen::Matrix4d emTcw_nid_last, emTwc_nid_last;
    cv::Mat mTcw_last = lFrame->mTcw_nid_.clone();
    
    if(mTcw_last.empty())
        return false;
    
    emTcw_nid_last<<mTcw_last.ptr<float>(0)[0], mTcw_last.ptr<float>(0)[1], mTcw_last.ptr<float>(0)[2], mTcw_last.ptr<float>(0)[3],
                    mTcw_last.ptr<float>(1)[0], mTcw_last.ptr<float>(1)[1], mTcw_last.ptr<float>(1)[2], mTcw_last.ptr<float>(1)[3],
                    mTcw_last.ptr<float>(2)[0], mTcw_last.ptr<float>(2)[1], mTcw_last.ptr<float>(2)[2], mTcw_last.ptr<float>(2)[3],
                    mTcw_last.ptr<float>(3)[0], mTcw_last.ptr<float>(3)[1], mTcw_last.ptr<float>(3)[2], mTcw_last.ptr<float>(3)[3];
    
    emTwc_nid_last = emTcw_nid_last.inverse();

    //assign value to a certain cell number
    //std::cout<<"size of points in first level "<<lFrame->hg_keys_pyr_[level].size()<<std::endl;
    for(int i = 0; i < lFrame->hg_keys_pyr_[level].size(); i++){
        int cell_num = lFrame->hg_keys_pyr_[level][i].class_id;
        if(cell_counter[cell_num]>one_nid_edge_size_-1 || cell_num == -1){
            continue;
        }

        float z_p = lFrame->keys_depth_pyr_[level][i];
        int pixel_u = lFrame->hg_keys_pyr_[level][i].pt.x;
        int pixel_v = lFrame->hg_keys_pyr_[level][i].pt.y;

        float x_p = z_p * (pixel_u - lFrame->cx_pyr_[level]) / lFrame->fx_pyr_[level];
        float y_p = z_p * (pixel_v - lFrame->cy_pyr_[level]) / lFrame->fy_pyr_[level];

        Eigen::Vector3d lp3d = (emTwc_nid_last*Eigen::Vector4d(x_p,y_p,z_p,1)).head(3);

        hg_keys_3d_last[cell_num].push_back(lp3d);
        //hg_intensity_last[cell_num].push_back(lFrame->image_pyr_[level].at<float>(lFrame->hg_keys_pyr_[level][i].pt.y, lFrame->hg_keys_pyr_[level][i].pt.x));

        hg_intensity_last(cell_counter[cell_num],cell_num) = lFrame->image_pyr_[level].at<float>(lFrame->hg_keys_pyr_[level][i].pt.y, lFrame->hg_keys_pyr_[level][i].pt.x);
        cell_counter[cell_num]++;

    }

    // for(int i = 0; i < hg_keys_3d_last.size(); i++){
    //     std::cout<<"size of hg intensity points "<<hg_keys_3d_last[i].size()<<std::endl;
    // }

    // std::cout<<"exit at optimizationNID "<<std::endl;
    // exit(0);
    
    //Project2LastFrame(pFrame, lFrame, Velocity, hg_keys_last, hg_intensity_last, hg_keys_3d, hg_keys_depth_current);

    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_300::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolver_6_300::PoseMatrixType>();

    g2o::BlockSolver_6_300 * solver_ptr = new g2o::BlockSolver_6_300(linearSolver);
    //solver_ptr->setLambda(50.0);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    solver->setUserLambdaInit(10000.0);//decrease the step
    //std::cout<<"current lamda "<<solver->currentLambda()<<std::endl;
    optimizer.setAlgorithm(solver);

    optimizer.setVerbose(true);

    int nInitialCorrespondences=0;

    //test use identity matrix to initialize
    //cv::Mat ini = cv::Mat::eye(4,4,CV_32F);

    // Set Frame vertex
    g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
    vSE3->setEstimate(Converter::toSE3Quat(T_cw));//Velocity.inv()
    vSE3->setId(0);
    vSE3->setFixed(false);
    optimizer.addVertex(vSE3);

    // Set MapPoint vertices
    const int N = lFrame->hg_keys_pyr_[level].size();

    vector<g2o::EdgeSE3ProjectIntensityOnlyPoseNID*> vpEdges;
    vector<size_t> vnIndexEdge;
    vector<bool> inlier;
    inlier.reserve(cell_size);
    vpEdges.reserve(cell_size);
    vnIndexEdge.reserve(cell_size);

    const float deltaMono = sqrt(3.1);//5.991

    {
    int counter = 0;
    for(int i=0; i<cell_size; i++)
    {
        if(cell_counter[i] < one_nid_edge_size_)
            continue;
                
        nInitialCorrespondences++;
       
        g2o::EdgeSE3ProjectIntensityOnlyPoseNID* e = new g2o::EdgeSE3ProjectIntensityOnlyPoseNID();
        e->image_current_ = pFrame->image_pyr_[level];
        e->x_world_set_ = hg_keys_3d_last[i];//hg_keys_3d_last[i].data();
        
        e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));

        // Eigen::MatrixXd test = hg_intensity_last.col(i);
        // counter = 0;
        // for(int t = 0; t < lFrame->hg_keys_pyr_[level].size(); t++){
        //     if(lFrame->hg_keys_pyr_[level][t].class_id == i){
        //         float z_p = lFrame->keys_depth_pyr_[level][t];
        //         int pixel_u = lFrame->hg_keys_pyr_[level][t].pt.x;
        //         int pixel_v = lFrame->hg_keys_pyr_[level][t].pt.y;

        //         float x_p = z_p * (pixel_u - lFrame->cx_pyr_[level]) / lFrame->fx_pyr_[level];
        //         float y_p = z_p * (pixel_v - lFrame->cy_pyr_[level]) / lFrame->fy_pyr_[level];

        //         Eigen::Vector3d lp3d = (emTwc_nid_last*Eigen::Vector4d(x_p,y_p,z_p,1)).head(3);

        //         std::cout<<"intensity for this point is "<<test(counter, 0)<<", it shoul also be "<<lFrame->image_pyr_[level].at<float>(lFrame->hg_keys_pyr_[level][t].pt.y, lFrame->hg_keys_pyr_[level][t].pt.x)<<std::endl;
        //         std::cout<<"this pixel's 3d world point \n"<<e->x_world_set_[counter]<<" \n it should also be \n"<<lp3d<<std::endl;
        //         counter++;
        //         if(counter > 10)
        //             break;
        //     }
        // }

        e->setMeasurement(hg_intensity_last.col(i));
        //std::cout<<"what happened here??? 1.... "<<one_nid_edge_size_<<std::endl;
        Eigen::MatrixXd info = Eigen::MatrixXd::Identity(1,1);
        e->setInformation(info);
        //std::cout<<"what happened here??? 2 "<<std::endl;
        e->set_bspline_relates(bs_,bs_der_,bs_degree_,bin_num_);
        //std::cout<<"PoseOptimizationNID loop............. counter "<<counter++<<std::endl;
        g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
        e->setRobustKernel(rk);
        rk->setDelta(deltaMono);

        e->fx_ = pFrame->fx_pyr_[level];
        e->fy_ = pFrame->fy_pyr_[level];
        e->cx_ = pFrame->cx_pyr_[level];
        e->cy_ = pFrame->cy_pyr_[level];

        e->computeHref();

        optimizer.addEdge(e);
        vpEdges.push_back(e);
        vnIndexEdge.push_back(i);
        inlier.push_back(true);
    }
    }

    std::cout<<"initial correspondence "<<nInitialCorrespondences<<std::endl;

    if(nInitialCorrespondences<1)
        return 0;

    // We perform 4 optimizations, after each optimization we classify observation as inlier/outlier
    // At the next optimization, outliers are not included, but at the end they can be classified as inliers again.
    const float chi2threshold[4]={10, 10, 10, 10};
    const int its[4]={10,10,10,10};    

    int nGood=0;
    for(size_t it=0; it<4; it++)
    {
        //std::cout<<"enter optimization ............. "<<it<<std::endl;
        vSE3->setEstimate(Converter::toSE3Quat(T_cw));
        
        optimizer.initializeOptimization(0);
        optimizer.optimize(its[it]);

        //std::cout<<"finish one iteration optimization............. "<<it<<std::endl;

        nGood=0;
        for(size_t i=0, iend=vpEdges.size(); i<iend; i++)
        {
            g2o::EdgeSE3ProjectIntensityOnlyPoseNID* e = vpEdges[i];

            if(inlier[i])
            {
                e->computeError();
            }
            else{
                continue;
            }

            const float chi2 = e->chi2();
            //std::cout<<"the chi 2 of current edge is "<<chi2<<std::endl;

            if(chi2>chi2threshold[it])
            {                
                e->setLevel(1);
                inlier[i] = false;
            }
            else
            {
                e->setLevel(0);
                nGood++;
            }
        }

        if(optimizer.edges().size()<1)
            break;
    }
    //std::cout<<"exit at optimization whole image "<<std::endl;
    //exit(0); 
    //std::cout<<"finish all optimization once"<<std::endl;
    //std::cout<<"\n"<<std::endl;   

    // // Recover optimized pose and return number of inliers
    g2o::VertexSE3Expmap* vSE3_recov = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(0));
    g2o::SE3Quat SE3quat_recov = vSE3_recov->estimate();
    cv::Mat pose = Converter::toCvMat(SE3quat_recov);

    //std::cout<<"optimized relative pose "<<pose.inv()<<std::endl;
    //std::cout<<"last frame nid pose "<<lFrame->mTcw_nid_<<std::endl;

    if(!lFrame->mTcw_nid_.empty()){
        pFrame->SetPoseNID(pose);
    }

    return 0;//nInitialCorrespondences-nBad;
}

/*Use whole image instead of picking some high gradient points*/
int Optimizer::PoseOptimizationNIDWholeImage(Frame *pFrame, Frame *lFrame, const cv::Mat& T_cw){
    int level = 0;
    int cell_size = pFrame->cell_  * pFrame->cell_;

    Eigen::Matrix4d emTcw_nid_last, emTwc_nid_last;
    cv::Mat mTcw_last = lFrame->mTcw_nid_.clone();

    //std::cout<<"timestamp of current frame and last frame "<<std::setprecision(16)<<pFrame->mTimeStamp<<","<<lFrame->mTimeStamp<<std::endl;
    
    if(mTcw_last.empty())
        return false;
    
    emTcw_nid_last<<mTcw_last.ptr<float>(0)[0], mTcw_last.ptr<float>(0)[1], mTcw_last.ptr<float>(0)[2], mTcw_last.ptr<float>(0)[3],
                    mTcw_last.ptr<float>(1)[0], mTcw_last.ptr<float>(1)[1], mTcw_last.ptr<float>(1)[2], mTcw_last.ptr<float>(1)[3],
                    mTcw_last.ptr<float>(2)[0], mTcw_last.ptr<float>(2)[1], mTcw_last.ptr<float>(2)[2], mTcw_last.ptr<float>(2)[3],
                    mTcw_last.ptr<float>(3)[0], mTcw_last.ptr<float>(3)[1], mTcw_last.ptr<float>(3)[2], mTcw_last.ptr<float>(3)[3];
    
    emTwc_nid_last = emTcw_nid_last.inverse();

    std::cout<<"current pose est \n"<<T_cw<<std::endl;
    // std::cout<<"last pose est \n"<<mTcw_last<<std::endl;
    // std::cout<<"last pose inverse \n"<<emTwc_nid_last<<std::endl;
    // std::cout<<"fx,fy,cx,cz "<<pFrame->fx_pyr_[level]<<","<<pFrame->fy_pyr_[level]<<","<<pFrame->cx_pyr_[level]<<","<<pFrame->cy_pyr_[level]<<std::endl;

    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_X::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolver_6_X::PoseMatrixType>();

    g2o::BlockSolver_6_X * solver_ptr = new g2o::BlockSolver_6_X(linearSolver);
    //solver_ptr->setLambda(50.0);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    //solver->setUserLambdaInit(100.0);//decrease the step
    //std::cout<<"current lamda "<<solver->currentLambda()<<std::endl;
    optimizer.setAlgorithm(solver);

    optimizer.setVerbose(true);

    int nInitialCorrespondences=0;

    //test use identity matrix to initialize
    cv::Mat ini = cv::Mat::eye(4,4,CV_32F);

    // Set Frame vertex
    g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
    vSE3->setEstimate(Converter::toSE3Quat(T_cw));//T_cw
    vSE3->setId(0);
    vSE3->setFixed(false);
    optimizer.addVertex(vSE3);

    vector<g2o::EdgeSE3ProjectIntensityOnlyPoseNID*> vpEdges;
    vector<size_t> vnIndexEdge;
    vector<bool> inlier;
    inlier.reserve(cell_size);
    vpEdges.reserve(cell_size);
    vnIndexEdge.reserve(cell_size);

    const float deltaMono = sqrt(0.8);

    {
        int counter = 0;
        for(int i=0; i<pFrame->cell_; i++){

            for(int j = 0; j<pFrame->cell_; j++)
            {
                g2o::EdgeSE3ProjectIntensityOnlyPoseNID* e = new g2o::EdgeSE3ProjectIntensityOnlyPoseNID();
                nInitialCorrespondences++;
                e->image_current_ = pFrame->image_pyr_[level];

                e->fx_ = pFrame->fx_pyr_[level];
                e->fy_ = pFrame->fy_pyr_[level];
                e->cx_ = pFrame->cx_pyr_[level];
                e->cy_ = pFrame->cy_pyr_[level];
                
                std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > point_3d_last;
                Eigen::Matrix<double, Eigen::Dynamic, 1> intensity_last;
                std::vector<cv::KeyPoint> point_2d_last;

                //Get3dPointAndIntensity(pFrame->cell_, i,j,lFrame->image_pyr_[level], lFrame->depth_pyr_[level], point_3d_last,intensity_last, emTwc_nid_last, pFrame->fx_pyr_[level], pFrame->fy_pyr_[level], pFrame->cx_pyr_[level], pFrame->cy_pyr_[level]);

                Get3dPointAndIntensity(pFrame->cell_, i,j,lFrame->image_pyr_[level], lFrame->depth_pyr_[level], point_3d_last,intensity_last, point_2d_last,emTwc_nid_last, pFrame->fx_pyr_[level], pFrame->fy_pyr_[level], pFrame->cx_pyr_[level], pFrame->cy_pyr_[level]);

                if(point_3d_last.size()<500)
                    continue;

                e->x_world_set_ = point_3d_last;
                e->hg_points_ = point_2d_last;
                
                e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));

                e->setMeasurement(intensity_last);

                Eigen::MatrixXd info = Eigen::MatrixXd::Identity(1,1);
                e->setInformation(info);

                e->set_bspline_relates(bs_,bs_der_,bs_degree_,bin_num_);

                g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                e->setRobustKernel(rk);
                rk->setDelta(deltaMono);

                e->computeHref();

                //std::cout<<"point 3d last size "<<point_3d_last.size()<<std::endl;

                optimizer.addEdge(e);
                vpEdges.push_back(e);
                vnIndexEdge.push_back(i);
                inlier.push_back(true);
            }

        }
    }

    //std::cout<<"size of edges "<<vpEdges.size()<<std::endl;

    const float chi2threshold[4]={0.8, 0.8, 0.99, 0.99};
    const int its[4]={10,10,10,10};    

    int nGood=0;
    for(size_t it=0; it<2; it++)
    {
        //std::cout<<"enter optimization ............. "<<it<<std::endl;
        vSE3->setEstimate(Converter::toSE3Quat(T_cw));
        
        optimizer.initializeOptimization(0);
        optimizer.optimize(its[it]);

        nGood=0;
        for(size_t i=0, iend=vpEdges.size(); i<iend; i++)
        {
            g2o::EdgeSE3ProjectIntensityOnlyPoseNID* e = vpEdges[i];

            if(inlier[i])
            {
                e->computeError();
            }
            else{
                continue;
            }

            const float chi2 = e->chi2();
            //std::cout<<"the chi 2 of current edge is "<<chi2<<std::endl;

            if(chi2>chi2threshold[it])
            {                
                e->setLevel(1);
                inlier[i] = false;
            }
            else
            {
                e->setLevel(0);
                nGood++;
            }
        }

        if(optimizer.edges().size()<1){
            std::cout<<"no enough edges"<<std::endl;
            break;
        }
    } 

    // Recover optimized pose and return number of inliers
    g2o::VertexSE3Expmap* vSE3_recov = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(0));
    g2o::SE3Quat SE3quat_recov = vSE3_recov->estimate();
    cv::Mat pose = Converter::toCvMat(SE3quat_recov);

    std::cout<<"pose optimized \n"<<pose<<std::endl;

    //std::cout<<"debug, test in pose optimization NID"<<std::endl;
    //exit(0);

    if(!lFrame->mTcw_nid_.empty()){
        pFrame->SetPoseNID(pose);
    }

    // std::cout<<"exit optimizer PoseOptimizationNIDWholeImage"<<std::endl;
    // exit(0);

    return 0;

}

void Optimizer::Get3dPointAndIntensity(int cell_size, int cell_row_id, int cell_col_id, const cv::Mat& image, const cv::Mat& depth, std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> >& points_3d, Eigen::VectorXd& intensity, const Eigen::Matrix4d& T_wc, float fx, float fy, float cx, float cy){
    int row_block = image.rows / cell_size;
    int col_block = image.cols / cell_size;
    int row_start = row_block * cell_row_id;
    int col_start = col_block * cell_col_id;
    int row_end = row_block * (cell_row_id + 1);
    int col_end = col_block * (cell_col_id + 1);
    int counter = 0;
    //std::cout<<"row start point and end point,  "<<row_start<<", "<<row_end<<", column start point and end point "<<col_start<<", "<<col_end<<std::endl;
    std::vector<double> intensity_v;
    for(int i = row_start; i < row_end; i++)
        for(int j = col_start; j < col_end; j++){
            float z_p = depth.ptr<float>(i)[j];
            counter++;
            if(z_p <0.01 || z_p > 100)
                continue;

            float x_p = z_p * (j - cx) / fx;
            float y_p = z_p * (i - cy) / fy;

            Eigen::Vector3d p_world = (T_wc*Eigen::Vector4d(x_p,y_p,z_p,1)).head(3);

            points_3d.push_back(p_world);
            intensity_v.push_back(image.ptr<float>(i)[j]);
            
        }
    
    intensity = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(intensity_v.data(), intensity_v.size());

    // for(int i = 0; i<counter; i++){
    //     std::cout<<"intensity "<<intensity(i)<<", "<<intensity_v[i]<<std::endl;
    // }

    //std::cout<<"counter is "<<counter<<", vector size is "<<intensity_v.size()<<std::endl;

    return;
};

void Optimizer::Get3dPointAndIntensity(int cell_size, int cell_row_id, int cell_col_id, const cv::Mat& image, const cv::Mat&depth, std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> >& points_3d, Eigen::Matrix<double, Eigen::Dynamic, 1>& intensity, std::vector<cv::KeyPoint>& points_2d, const Eigen::Matrix4d& T_wc, float fx, float fy, float cx, float cy){
    int row_block = image.rows / cell_size;
    int col_block = image.cols / cell_size;
    int row_start = row_block * cell_row_id;
    int col_start = col_block * cell_col_id;
    int row_end = row_block * (cell_row_id + 1);
    int col_end = col_block * (cell_col_id + 1);
    int counter = 0;
    //std::cout<<"row start point and end point,  "<<row_start<<", "<<row_end<<", column start point and end point "<<col_start<<", "<<col_end<<std::endl;
    std::vector<double> intensity_v;
    for(int i = row_start; i < row_end; i++)
        for(int j = col_start; j < col_end; j++){
            float z_p = depth.ptr<float>(i)[j];
            counter++;
            if(z_p <0.01 || z_p > 100)
                continue;

            float x_p = z_p * (j - cx) / fx;
            float y_p = z_p * (i - cy) / fy;

            Eigen::Vector3d p_world = (T_wc*Eigen::Vector4d(x_p,y_p,z_p,1)).head(3);

            // if(counter%5 == 0)
            //     std::cout<<"original pixel x,y,depth "<<j<<","<<i<<","<<z_p<<std::endl;

            cv::KeyPoint pt2d;
            pt2d.pt.x = j;
            pt2d.pt.y = i;

            points_2d.push_back(pt2d); 
            points_3d.push_back(p_world);
            intensity_v.push_back(image.ptr<float>(i)[j]);
            
        }
    //std::cout<<"size of vector intensity, 3d point "<<intensity_v.size()<<","<<points_3d.size()<<std::endl;
    intensity = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(intensity_v.data(), intensity_v.size());
    //std::cout<<"size of eigen intensity"<<intensity.rows()<<std::endl;
    return;    
}


// int Optimizer::PoseOptimizationIntensity(Frame *pFrame, Frame *lFrame, const cv::Mat& Velocity){

//     //std::cout<<"enter pose optimization intensity"<<std::endl;

//     g2o::SparseOptimizer optimizer;
//     g2o::BlockSolver_6_1::LinearSolverType * linearSolver;

//     linearSolver = new g2o::LinearSolverDense<g2o::BlockSolver_6_1::PoseMatrixType>();

//     g2o::BlockSolver_6_1 * solver_ptr = new g2o::BlockSolver_6_1(linearSolver);
//     //solver_ptr->setLambda(50.0);

//     g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
//     optimizer.setAlgorithm(solver);

//     //optimizer.setVerbose(true);

//     int nInitialCorrespondences=0;

//     // Set Frame vertex
//     g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();

//     //test use identity matrix to initialize
//     cv::Mat ini = cv::Mat::eye(4,4,CV_32F);
    
//     //std::cout<<"before optimization the matrix is \n "<<Velocity.inv()<<std::endl;
//     vSE3->setEstimate(Converter::toSE3Quat(Velocity.inv()));
//     vSE3->setId(0);
//     vSE3->setFixed(false);
//     optimizer.addVertex(vSE3);

//     int counter = 0;
//     int pyr_num = 6;
//     int N = 0;
//     // Set MapPoint vertices
//     for(int i = 0; i<pyr_num; i++){
//         N += pFrame->hg_keys_pyr_[i].size();
//     }
    
//     vector<g2o::EdgeSE3ProjectIntensityOnlyPose*> vpEdges;
//     vector<size_t> vnIndexEdge;
//     vector<bool> inlier(N, true);
//     vpEdges.reserve(N);
//     vnIndexEdge.reserve(N);

//     const float delta = sqrt(50);//5.991

//     // for(int i = 0; i<6; i++)
//     //     std::cout<<"camera parameters "<<pFrame->fx_pyr_[i]<<","<<pFrame->fy_pyr_[i]<<","<<pFrame->cx_pyr_[i]<<","<<pFrame->cy_pyr_[i]<<std::endl;

//     {
//     // unique_lock<mutex> lock(MapPoint::mGlobalMutex);

//     std::vector<int> start_index(pyr_num, 0);
//     for(int level = 1; level < pyr_num; level++){
//         start_index[level] += start_index[level-1] + pFrame->hg_keys_pyr_[level-1].size();  
//     }
//     for(int level = 0; level < pyr_num; level++){
//         //std::cout<<"the hg keys number of level "<<level<< " is "<<pFrame->hg_keys_pyr_[level].size()<<std::endl;
//         for(int i=0; i<pFrame->hg_keys_pyr_[level].size(); i++)
//         {
//             nInitialCorrespondences++;

//             g2o::EdgeSE3ProjectIntensityOnlyPose* e = new g2o::EdgeSE3ProjectIntensityOnlyPose();
//             e->image_last_ = lFrame->image_pyr_[level];
//             e->hg_point_current_ = pFrame->hg_keys_pyr_[level][i];
//             e->hg_depth_current_ = pFrame->keys_depth_pyr_[level][i];

//             e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
//             Eigen::Matrix<double, 1,1> current_intensity;
//             cv::Mat image_level = pFrame->image_pyr_[level];
//             int position_x = pFrame->hg_keys_pyr_[level][i].pt.x;
//             int position_y = pFrame->hg_keys_pyr_[level][i].pt.y;
//             current_intensity(0,0) = image_level.ptr<float>(position_y)[position_x];
//             e->setMeasurement(current_intensity);

//             //lower level has less points so more weights
//             float weight = pow(1.5,level);//2 * level + 1;
//             Eigen::MatrixXd info = Eigen::MatrixXd::Identity(1,1);
//             e->setInformation(weight*info);

//             //std::cout<<"PoseOptimizationNID loop............. counter "<<counter++<<std::endl;
//             g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
//             e->setRobustKernel(rk);
//             rk->setDelta(delta);

//             e->fx_ = pFrame->fx_pyr_[level];
//             e->fy_ = pFrame->fy_pyr_[level];
//             e->cx_ = pFrame->cx_pyr_[level];
//             e->cy_ = pFrame->cy_pyr_[level];

//             optimizer.addEdge(e);
//             vpEdges.push_back(e);
//             //vnIndexEdge.push_back(i);
//             vnIndexEdge.push_back(start_index[level] + i);
//         }
//     }
    
//     }

//     //std::cout<<"raw edge number "<<vpEdges.size()<<std::endl;

//     // double all_chi2 = 0;
//     // for(size_t i=0, iend=vpEdges.size(); i<iend; i++){
//     //     g2o::EdgeSE3ProjectIntensityOnlyPose* e = vpEdges[i];
//     //     e->computeError();
//     //     const double chi2 = e->chi2();
//     //     all_chi2 += chi2;
//     // }
//     // all_frame_counter_++;
//     // all_accu_error_ += all_chi2;
//     // std::cout<<"all frame counter and all accumulate error and average"<<std::setprecision(10)<<all_frame_counter_<<", "<<all_accu_error_<<","<<all_accu_error_/all_frame_counter_<<std::endl;
//     // std::cout<<"initial chi2 error is "<<all_chi2<<std::endl;

//     if(nInitialCorrespondences<3)
//         return 0;

//     // We perform 4 optimizations, after each optimization we classify observation as inlier/outlier
//     // At the next optimization, outliers are not included, but at the end they can be classified as inliers again.
//     const float chi2threshold[4]={50, 50, 50, 50};
//     const int its[4]={10,10,10,10};    

//     for(size_t it=0; it<4; it++)
//     {
//         vSE3->setEstimate(Converter::toSE3Quat(Velocity.inv()));

//         //std::cout<<"enter optimization............................... "<<it<<std::endl;

//         optimizer.initializeOptimization(0);
//         optimizer.optimize(its[it]);

//         //std::cout<<"finish optimization........................... "<<it<<std::endl;
//         //std::cout<<std::endl;

//         for(size_t i=0, iend=vpEdges.size(); i<iend; i++)
//         {
//             g2o::EdgeSE3ProjectIntensityOnlyPose* e = vpEdges[i];

//             if(inlier[i])
//                 e->computeError();
//             else
//                 continue;

//             const float chi2 = e->chi2();

//             if(chi2>chi2threshold[it])
//             {                
//                 e->setLevel(1);
//             }
//             else
//             {
//                 e->setLevel(0);
//                 inlier[i] = false;
//             }

//             if(it==2)
//                 e->setRobustKernel(0);

//         }

//         if(optimizer.edges().size()<10)
//             break;
//     }

//     // all_frame_counter_++;
//     // std::cout<<"frame number is "<<all_frame_counter_<<", current all average edges "<<all_available_edge_/all_frame_counter_<<std::endl;

//     // if(all_available_edge_>350){
//     //     all_frame_counter_++;
//     //     std::cout<<"there are "<<all_frame_counter_<<" have edge number bigger than "<<350<<std::endl;
//     // }
//     // all_available_edge_ = 0;

//     //Recover optimized pose and return number of inliers
//     g2o::VertexSE3Expmap* vSE3_recov = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(0));
//     g2o::SE3Quat SE3quat_recov = vSE3_recov->estimate();
//     cv::Mat pose = Converter::toCvMat(SE3quat_recov);

//     if(!lFrame->mTcw_pho_.empty()){
//         pFrame->SetPosePho(pose.inv() * lFrame->mTcw_pho_);
//     }

//     //std::cout<<"after optimization the matrix is \n "<<pose<<std::endl;
//     std::cout<<" \n \n"<<std::endl;

//     return 0;
// }

// float Optimizer::PoseOptimizationIntensityCoarse(Frame *pFrame, Frame *lFrame, cv::Mat& Velocity){
//     //std::cout<<"enter pose optimization intensity"<<std::endl;

//     g2o::SparseOptimizer optimizer;
//     g2o::BlockSolver_6_1::LinearSolverType * linearSolver;

//     linearSolver = new g2o::LinearSolverDense<g2o::BlockSolver_6_1::PoseMatrixType>();

//     g2o::BlockSolver_6_1 * solver_ptr = new g2o::BlockSolver_6_1(linearSolver);
//     //solver_ptr->setLambda(50.0);

//     g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
//     optimizer.setAlgorithm(solver);

//     optimizer.setVerbose(true);

//     int nInitialCorrespondences=0;

//     // Set Frame vertex
//     g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();

//     //test use identity matrix to initialize
//     cv::Mat ini = cv::Mat::eye(4,4,CV_32F);
    
//     //std::cout<<"before optimization the matrix is \n "<<Velocity.inv()<<std::endl;
//     vSE3->setEstimate(Converter::toSE3Quat(Velocity.inv()));
//     vSE3->setId(0);
//     vSE3->setFixed(false);
//     optimizer.addVertex(vSE3);

//     int level = 0;//use 1 for less points
//     // Set MapPoint vertices
//     const int N = pFrame->hg_keys_pyr_[level].size();

//     std::cout<<"number of high gradient points in optimization "<<N<<std::endl;

//     vector<g2o::EdgeSE3ProjectIntensityOnlyPose*> vpEdges;
//     vector<size_t> vnIndexEdge;
//     vector<bool> inlier(N, true);
//     vpEdges.reserve(N);
//     vnIndexEdge.reserve(N);

//     const float delta = sqrt(50);//5.991

//     // for(int i = 0; i<6; i++)
//     //     std::cout<<"camera parameters "<<pFrame->fx_pyr_[i]<<","<<pFrame->fy_pyr_[i]<<","<<pFrame->cx_pyr_[i]<<","<<pFrame->cy_pyr_[i]<<std::endl;

//     {
//     // unique_lock<mutex> lock(MapPoint::mGlobalMutex);
//     int counter = 0;
//     int pyr_num = 6;
//     std::vector<int> start_index(pyr_num, 0);

//     //std::cout<<"the hg keys number of level "<<level<< " is "<<pFrame->hg_keys_pyr_[level].size()<<std::endl;
//     for(int i=0; i<pFrame->hg_keys_pyr_[level].size(); i++)
//     {
//         nInitialCorrespondences++;

//         g2o::EdgeSE3ProjectIntensityOnlyPose* e = new g2o::EdgeSE3ProjectIntensityOnlyPose();
//         e->image_last_ = lFrame->image_pyr_[level];
//         e->hg_point_current_ = pFrame->hg_keys_pyr_[level][i];
//         e->hg_depth_current_ = pFrame->keys_depth_pyr_[level][i];

//         e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
//         Eigen::Matrix<double, 1,1> current_intensity;
//         cv::Mat image_level = pFrame->image_pyr_[level];
//         int position_x = pFrame->hg_keys_pyr_[level][i].pt.x;
//         int position_y = pFrame->hg_keys_pyr_[level][i].pt.y;
//         current_intensity(0,0) = image_level.ptr<float>(position_y)[position_x];
//         e->setMeasurement(current_intensity);

//         //lower level has less points so more weights
//         float weight = pow(1.5,level);//2 * level + 1;
//         Eigen::MatrixXd info = Eigen::MatrixXd::Identity(1,1);
//         e->setInformation(weight*info);

//         //std::cout<<"PoseOptimizationNID loop............. counter "<<counter++<<std::endl;
//         g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
//         e->setRobustKernel(rk);
//         rk->setDelta(delta);

//         e->fx_ = pFrame->fx_pyr_[level];
//         e->fy_ = pFrame->fy_pyr_[level];
//         e->cx_ = pFrame->cx_pyr_[level];
//         e->cy_ = pFrame->cy_pyr_[level];

//         optimizer.addEdge(e);
//         vpEdges.push_back(e);
//         vnIndexEdge.push_back(i);
//     }
    
    
//     }

//     //std::cout<<"raw edge number "<<vpEdges.size()<<std::endl;

//     // double all_chi2 = 0;
//     // for(size_t i=0, iend=vpEdges.size(); i<iend; i++){
//     //     g2o::EdgeSE3ProjectIntensityOnlyPose* e = vpEdges[i];
//     //     e->computeError();
//     //     const double chi2 = e->chi2();
//     //     all_chi2 += chi2;
//     // }
//     // all_frame_counter_++;
//     // all_accu_error_ += all_chi2;
//     // std::cout<<"all frame counter and all accumulate error and average"<<std::setprecision(10)<<all_frame_counter_<<", "<<all_accu_error_<<","<<all_accu_error_/all_frame_counter_<<std::endl;
//     // std::cout<<"initial chi2 error is "<<all_chi2<<std::endl;

//     if(nInitialCorrespondences<3)
//         return 0;

//     // We perform 4 optimizations, after each optimization we classify observation as inlier/outlier
//     // At the next optimization, outliers are not included, but at the end they can be classified as inliers again.
//     const float chi2threshold[2]={50, 50};
//     const int its[2]={10,10};    

//     int nGood=0;
//     for(size_t it=0; it<2; it++)
//     {
//         vSE3->setEstimate(Converter::toSE3Quat(Velocity.inv()));

//         //std::cout<<"enter optimization............................... "<<it<<std::endl;

//         optimizer.initializeOptimization(0);
//         optimizer.optimize(its[it]);

//         //std::cout<<"finish optimization........................... "<<it<<std::endl;
//         //std::cout<<std::endl;

//         nGood=0;
//         for(size_t i=0, iend=vpEdges.size(); i<iend; i++)
//         {
//             g2o::EdgeSE3ProjectIntensityOnlyPose* e = vpEdges[i];

//             if(inlier[i])
//             {
//                 e->computeError();
//             }
//             else{
//                 continue;
//             }

//             const float chi2 = e->chi2();

//             if(chi2>chi2threshold[it])
//             {                
//                 e->setLevel(1);
//                 inlier[i] = false;
//             }
//             else
//             {
//                 e->setLevel(0);
//                 if(it==1)
//                     nGood++;
//             }

//         }

//         //std::cout<<"bad edges "<<nBad<<std::endl;

//         if(optimizer.edges().size()<10)
//             break;
//     }

//     // all_frame_counter_++;
//     // std::cout<<"frame number is "<<all_frame_counter_<<", current all average edges "<<all_available_edge_/all_frame_counter_<<std::endl;

//     // if(all_available_edge_>350){
//     //     all_frame_counter_++;
//     //     std::cout<<"there are "<<all_frame_counter_<<" have edge number bigger than "<<350<<std::endl;
//     // }
//     // all_available_edge_ = 0;

//     //Recover optimized pose and return number of inliers
//     g2o::VertexSE3Expmap* vSE3_recov = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(0));
//     g2o::SE3Quat SE3quat_recov = vSE3_recov->estimate();
//     cv::Mat pose = Converter::toCvMat(SE3quat_recov);

//     //std::cout<<" \n \n"<<std::endl;

//     Velocity = pose.clone().inv();

//     return ((float)nGood)/vpEdges.size();

//     //std::cout<<"after optimization the matrix is \n "<<pose<<std::endl;
    
// };


bool Optimizer::PoseOptimizationDirectCoarse(Frame *pFrame, Frame *lFrame, cv::Mat T_cw){
    //std::cout<<"enter pose optimization intensity"<<std::endl;
    //cv::Mat T_cw = pFrame->mTcw_pho_.clone();

    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_1::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolver_6_1::PoseMatrixType>();

    g2o::BlockSolver_6_1 * solver_ptr = new g2o::BlockSolver_6_1(linearSolver);
    //solver_ptr->setLambda(50.0);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    //optimizer.setVerbose(true);

    int nInitialCorrespondences=0;

    // Set Frame vertex
    g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();

    //test use identity matrix to initialize
    cv::Mat ini = cv::Mat::eye(4,4,CV_32F);
    
    //std::cout<<"before optimization the matrix is \n "<<T_cw<<std::endl;
    vSE3->setEstimate(Converter::toSE3Quat(T_cw));
    vSE3->setId(0);
    vSE3->setFixed(false);
    optimizer.addVertex(vSE3);

    int level = 1;//use 1 for less points
    // Set MapPoint vertices
    const int N = lFrame->hg_keys_pyr_[level].size();

    //std::cout<<"number of high gradient points in optimization "<<N<<std::endl;

    vector<g2o::EdgeSE3ProjectDirectOnlyPose*> vpEdges;
    vector<size_t> vnIndexEdge;
    vector<bool> inlier(N, true);
    vpEdges.reserve(N);
    vnIndexEdge.reserve(N);

    Eigen::Matrix4d emTcw_pho_last, emTwc_pho_last;
    cv::Mat mTcw_last = lFrame->mTcw_pho_.clone();
    if(mTcw_last.empty())
        return false;
    emTcw_pho_last<<mTcw_last.ptr<float>(0)[0], mTcw_last.ptr<float>(0)[1], mTcw_last.ptr<float>(0)[2], mTcw_last.ptr<float>(0)[3],
                    mTcw_last.ptr<float>(1)[0], mTcw_last.ptr<float>(1)[1], mTcw_last.ptr<float>(1)[2], mTcw_last.ptr<float>(1)[3],
                    mTcw_last.ptr<float>(2)[0], mTcw_last.ptr<float>(2)[1], mTcw_last.ptr<float>(2)[2], mTcw_last.ptr<float>(2)[3],
                    mTcw_last.ptr<float>(3)[0], mTcw_last.ptr<float>(3)[1], mTcw_last.ptr<float>(3)[2], mTcw_last.ptr<float>(3)[3];
    
    emTwc_pho_last = emTcw_pho_last.inverse();
    const float delta = sqrt(100);//5.991

    // for(int i = 0; i<6; i++)
    //     std::cout<<"camera parameters "<<pFrame->fx_pyr_[i]<<","<<pFrame->fy_pyr_[i]<<","<<pFrame->cx_pyr_[i]<<","<<pFrame->cy_pyr_[i]<<std::endl;

    {
    // unique_lock<mutex> lock(MapPoint::mGlobalMutex);
    int counter = 0;
    int pyr_num = 6;

    //std::cout<<"the hg keys number of level "<<level<< " is "<<pFrame->hg_keys_pyr_[level].size()<<std::endl;
    for(int i=0; i<lFrame->hg_keys_pyr_[level].size(); i++)
    {
        nInitialCorrespondences++;

        g2o::EdgeSE3ProjectDirectOnlyPose* e = new g2o::EdgeSE3ProjectDirectOnlyPose();
        e->image_current_ = pFrame->image_pyr_[level];

        float z_p = lFrame->keys_depth_pyr_[level][i];
        int pixel_u = lFrame->hg_keys_pyr_[level][i].pt.x;
        int pixel_v = lFrame->hg_keys_pyr_[level][i].pt.y;

        float x_p = z_p * (pixel_u - lFrame->cx_pyr_[level]) / lFrame->fx_pyr_[level];
        float y_p = z_p * (pixel_v - lFrame->cy_pyr_[level]) / lFrame->fy_pyr_[level];

        e->x_world_ = (emTwc_pho_last*Eigen::Vector4d(x_p,y_p,z_p,1)).head(3);
         
        e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
        double last_intensity;
        cv::Mat image_level = lFrame->image_pyr_[level];
        last_intensity = image_level.ptr<float>(pixel_v)[pixel_u];
        e->setMeasurement(last_intensity);

        //lower level has less points so more weights
        float weight = 1;//pow(1.5,level);//2 * level + 1;
        Eigen::MatrixXd info = Eigen::MatrixXd::Identity(1,1);
        e->setInformation(weight*info);

        //std::cout<<"PoseOptimizationNID loop............. counter "<<counter++<<std::endl;
        g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
        e->setRobustKernel(rk);
        rk->setDelta(delta);

        e->fx_ = lFrame->fx_pyr_[level];
        e->fy_ = lFrame->fy_pyr_[level];
        e->cx_ = lFrame->cx_pyr_[level];
        e->cy_ = lFrame->cy_pyr_[level];

        optimizer.addEdge(e);
        vpEdges.push_back(e);
        vnIndexEdge.push_back(i);
    }
    
    
    }

    //std::cout<<"raw edge number "<<vpEdges.size()<<std::endl;

    // double all_chi2 = 0;
    // for(size_t i=0, iend=vpEdges.size(); i<iend; i++){
    //     g2o::EdgeSE3ProjectIntensityOnlyPose* e = vpEdges[i];
    //     e->computeError();
    //     const double chi2 = e->chi2();
    //     all_chi2 += chi2;
    // }
    // all_frame_counter_++;
    // all_accu_error_ += all_chi2;
    // std::cout<<"all frame counter and all accumulate error and average"<<std::setprecision(10)<<all_frame_counter_<<", "<<all_accu_error_<<","<<all_accu_error_/all_frame_counter_<<std::endl;
    // std::cout<<"initial chi2 error is "<<all_chi2<<std::endl;

    if(nInitialCorrespondences<3)
        return false;

    // We perform 4 optimizations, after each optimization we classify observation as inlier/outlier
    // At the next optimization, outliers are not included, but at the end they can be classified as inliers again.
    const float chi2threshold[2]={100, 100};
    const int its[2]={5,10};    

    int nGood=0;
    for(size_t it=0; it<1; it++)
    {
        vSE3->setEstimate(Converter::toSE3Quat(T_cw));

        //std::cout<<"enter optimization............................... "<<it<<std::endl;

        optimizer.initializeOptimization(0);
        int it_num = optimizer.optimize(its[it]);

        if(it_num <= 1)
            return false;

        //std::cout<<"finish optimization........................... "<<it<<std::endl;
        const double* chi2_his = optimizer.activeRobustChi2His();

        double decrease_rate = (chi2_his[0] - chi2_his[it_num-1])/chi2_his[0];

        if(decrease_rate < 0.1)
            return false;

        // nGood=0;
        // for(size_t i=0; i<vpEdges.size(); i++)
        // {
        //     g2o::EdgeSE3ProjectDirectOnlyPose* e = vpEdges[i];

        //     if(inlier[i])
        //     {
        //         e->computeError();
        //     }
        //     else{
        //         continue;
        //     }

        //     const float chi2 = e->chi2();

        //     if(chi2>chi2threshold[it])
        //     {                
        //         e->setLevel(1);
        //         inlier[i] = false;
        //     }
        //     else
        //     {
        //         e->setLevel(0);
        //         nGood++;
        //     }

        // }

        //std::cout<<"good edges "<<nGood<<std::endl;

        if(optimizer.edges().size()<10){
            return false;
        }
    }

    // all_frame_counter_++;
    // std::cout<<"frame number is "<<all_frame_counter_<<", current all average edges "<<all_available_edge_/all_frame_counter_<<std::endl;

    // if(all_available_edge_>350){
    //     all_frame_counter_++;
    //     std::cout<<"there are "<<all_frame_counter_<<" have edge number bigger than "<<350<<std::endl;
    // }
    // all_available_edge_ = 0;

    //Recover optimized pose and return number of inliers
    // g2o::VertexSE3Expmap* vSE3_recov = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(0));
    // g2o::SE3Quat SE3quat_recov = vSE3_recov->estimate();
    // cv::Mat pose = Converter::toCvMat(SE3quat_recov);

    //std::cout<<" \n \n"<<std::endl;
    //std::cout<<"after optimization the pose is \n "<<pose<<std::endl;
    
    //pFrame->SetPosePho(pose);

    return true;
    
};

bool Optimizer::PoseOptimizationDirect(Frame *pFrame, Frame *lFrame, cv::Mat T_cw){

    //cv::Mat T_cw = pFrame->mTcw_pho_.clone();

    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_1::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolver_6_1::PoseMatrixType>();

    g2o::BlockSolver_6_1 * solver_ptr = new g2o::BlockSolver_6_1(linearSolver);
    //solver_ptr->setLambda(50.0);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    optimizer.setVerbose(true);

    int nInitialCorrespondences=0;

    // Set Frame vertex
    g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();

    //test use identity matrix to initialize
    cv::Mat ini = cv::Mat::eye(4,4,CV_32F);
    
    //std::cout<<"before optimization the matrix is \n "<<T_cw<<std::endl;
    vSE3->setEstimate(Converter::toSE3Quat(T_cw));
    vSE3->setId(0);
    vSE3->setFixed(false);
    optimizer.addVertex(vSE3);

    // Set MapPoint vertices
    int N = 0;
    int pyr_num = 6;
    // Set MapPoint vertices
    for(int i = 0; i<pyr_num; i++){
        N += lFrame->hg_keys_pyr_[i].size();
    }

    std::vector<int> start_index(pyr_num, 0);
    for(int level = 1; level < pyr_num; level++){
        start_index[level] += start_index[level-1] + lFrame->hg_keys_pyr_[level-1].size();  
    }

    //std::cout<<"number of high gradient points in optimization "<<N<<std::endl;

    vector<g2o::EdgeSE3ProjectDirectOnlyPose*> vpEdges;
    vector<size_t> vnIndexEdge;
    vector<bool> inlier(N, true);
    vpEdges.reserve(N);
    vnIndexEdge.reserve(N);

    Eigen::Matrix4d emTcw_pho_last, emTwc_pho_last;
    cv::Mat mTcw_last = lFrame->mTcw_pho_.clone();
    if(mTcw_last.empty())
        return false;
    emTcw_pho_last<<mTcw_last.ptr<float>(0)[0], mTcw_last.ptr<float>(0)[1], mTcw_last.ptr<float>(0)[2], mTcw_last.ptr<float>(0)[3],
                    mTcw_last.ptr<float>(1)[0], mTcw_last.ptr<float>(1)[1], mTcw_last.ptr<float>(1)[2], mTcw_last.ptr<float>(1)[3],
                    mTcw_last.ptr<float>(2)[0], mTcw_last.ptr<float>(2)[1], mTcw_last.ptr<float>(2)[2], mTcw_last.ptr<float>(2)[3],
                    mTcw_last.ptr<float>(3)[0], mTcw_last.ptr<float>(3)[1], mTcw_last.ptr<float>(3)[2], mTcw_last.ptr<float>(3)[3];
    
    emTwc_pho_last = emTcw_pho_last.inverse();
    const float delta = sqrt(100);//5.991

    // for(int i = 0; i<6; i++)
    //     std::cout<<"camera parameters "<<pFrame->fx_pyr_[i]<<","<<pFrame->fy_pyr_[i]<<","<<pFrame->cx_pyr_[i]<<","<<pFrame->cy_pyr_[i]<<std::endl;

    {
    // unique_lock<mutex> lock(MapPoint::mGlobalMutex);
    int counter = 0;

    //std::cout<<"the hg keys number of level "<<level<< " is "<<pFrame->hg_keys_pyr_[level].size()<<std::endl;
    for(int level = 0; level < pyr_num; level++){
        for(int i=0; i<lFrame->hg_keys_pyr_[level].size(); i++)
        {
            nInitialCorrespondences++;

            g2o::EdgeSE3ProjectDirectOnlyPose* e = new g2o::EdgeSE3ProjectDirectOnlyPose();
            e->image_current_ = pFrame->image_pyr_[level];

            float z_p = lFrame->keys_depth_pyr_[level][i];
            int pixel_u = lFrame->hg_keys_pyr_[level][i].pt.x;
            int pixel_v = lFrame->hg_keys_pyr_[level][i].pt.y;

            float x_p = z_p * (pixel_u - lFrame->cx_pyr_[level]) / lFrame->fx_pyr_[level];
            float y_p = z_p * (pixel_v - lFrame->cy_pyr_[level]) / lFrame->fy_pyr_[level];

            e->x_world_ = (emTwc_pho_last*Eigen::Vector4d(x_p,y_p,z_p,1)).head(3);
            
            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
            double last_intensity;
            cv::Mat image_level = lFrame->image_pyr_[level];
            last_intensity = image_level.ptr<float>(pixel_v)[pixel_u];
            e->setMeasurement(last_intensity);

            //lower level has less points so more weights
            float weight = 1;//pow(1.5,level);//2 * level + 1;
            Eigen::MatrixXd info = Eigen::MatrixXd::Identity(1,1);
            e->setInformation(weight*info);

            //std::cout<<"PoseOptimizationNID loop............. counter "<<counter++<<std::endl;
            g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
            e->setRobustKernel(rk);
            rk->setDelta(delta);

            e->fx_ = lFrame->fx_pyr_[level];
            e->fy_ = lFrame->fy_pyr_[level];
            e->cx_ = lFrame->cx_pyr_[level];
            e->cy_ = lFrame->cy_pyr_[level];

            optimizer.addEdge(e);
            vpEdges.push_back(e);
            vnIndexEdge.push_back(start_index[level] + i);
        }
    }
    
    
    }

    // std::cout<<"raw edge number "<<vpEdges.size()<<std::endl;
    // std::cout<<"inlier index size "<<inlier.size()<<std::endl;

    // double all_chi2 = 0;
    // for(size_t i=0, iend=vpEdges.size(); i<iend; i++){
    //     g2o::EdgeSE3ProjectIntensityOnlyPose* e = vpEdges[i];
    //     e->computeError();
    //     const double chi2 = e->chi2();
    //     all_chi2 += chi2;
    // }
    // all_frame_counter_++;
    // all_accu_error_ += all_chi2;
    // std::cout<<"all frame counter and all accumulate error and average"<<std::setprecision(10)<<all_frame_counter_<<", "<<all_accu_error_<<","<<all_accu_error_/all_frame_counter_<<std::endl;
    // std::cout<<"initial chi2 error is "<<all_chi2<<std::endl;

    if(nInitialCorrespondences<3)
        return false;

    // We perform 4 optimizations, after each optimization we classify observation as inlier/outlier
    // At the next optimization, outliers are not included, but at the end they can be classified as inliers again.
    const float chi2threshold[4]={100, 100, 100, 100};
    const int its[4]={15,10,10,10};    

    int nGood=0;
    for(size_t it=0; it<2; it++)
    {
        vSE3->setEstimate(Converter::toSE3Quat(T_cw));

        //std::cout<<"enter optimization............................... "<<it<<std::endl;

        optimizer.initializeOptimization(0);
        optimizer.optimize(its[it]);

        //std::cout<<"finish optimization........................... "<<it<<std::endl;

        // nGood=0;
        for(size_t i=0; i<vpEdges.size(); i++)
        {
            g2o::EdgeSE3ProjectDirectOnlyPose* e = vpEdges[i];

            if(inlier[i])
            {
                e->computeError();
            }
            else{
                continue;
            }

            const float chi2 = e->chi2();

            if(chi2>chi2threshold[it])
            {                
                e->setLevel(1);
                inlier[i] = false;
            }
            else
            {
                e->setLevel(0);
                nGood++;
            }

        }

        //std::cout<<"good edges "<<nGood<<std::endl;

        if(optimizer.edges().size()<10)
            break;
    }

    // all_frame_counter_++;
    // std::cout<<"frame number is "<<all_frame_counter_<<", current all average edges "<<all_available_edge_/all_frame_counter_<<std::endl;

    //Recover optimized pose and return number of inliers
    g2o::VertexSE3Expmap* vSE3_recov = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(0));
    g2o::SE3Quat SE3quat_recov = vSE3_recov->estimate();
    cv::Mat pose = Converter::toCvMat(SE3quat_recov);

    // std::cout<<"after optimization the pose is \n "<<pose<<std::endl;
    // std::cout<<" \n \n"<<std::endl;
    
    pFrame->SetPosePho(pose);

    return true;
}

} //namespace ORB_SLAM
