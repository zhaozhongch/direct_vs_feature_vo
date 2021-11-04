// g2o - General Graph Optimization
// Copyright (C) 2011 H. Strasdat
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// Modified by Raúl Mur Artal (2014)
// Added EdgeSE3ProjectXYZ (project using focal_length in x,y directions)
// Modified by Raúl Mur Artal (2016)
// Added EdgeStereoSE3ProjectXYZ (project using focal_length in x,y directions)
// Added EdgeSE3ProjectXYZOnlyPose (unary edge to optimize only the camera pose)
// Added EdgeStereoSE3ProjectXYZOnlyPose (unary edge to optimize only the camera pose)

#ifndef G2O_SIX_DOF_TYPES_EXPMAP
#define G2O_SIX_DOF_TYPES_EXPMAP

#include "../core/base_vertex.h"
#include "../core/base_binary_edge.h"
#include "../core/base_unary_edge.h"
#include "se3_ops.h"
#include "se3quat.h"
#include "types_sba.h"
#include <Eigen/Geometry>
#include "opencv2/core.hpp"

namespace g2o {

namespace types_six_dof_expmap {
void init();
}

using namespace Eigen;

typedef Matrix<double, 6, 6> Matrix6d;

/**
 * \brief SE3 Vertex parameterized internally with a transformation matrix
 and externally with its exponential map
 */
class  VertexSE3Expmap : public BaseVertex<6, SE3Quat>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  VertexSE3Expmap();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  virtual void setToOriginImpl() {
    _estimate = SE3Quat();
  }

  virtual void oplusImpl(const double* update_)  {
    Eigen::Map<const Vector6d> update(update_);
    setEstimate(SE3Quat::exp(update)*estimate());
  }
};


class  EdgeSE3ProjectXYZ: public  BaseBinaryEdge<2, Vector2d, VertexSBAPointXYZ, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeSE3ProjectXYZ();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
    const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
    Vector2d obs(_measurement);
    _error = obs-cam_project(v1->estimate().map(v2->estimate()));
  }

  bool isDepthPositive() {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
    const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
    return (v1->estimate().map(v2->estimate()))(2)>0.0;
  }
    

  virtual void linearizeOplus();

  Vector2d cam_project(const Vector3d & trans_xyz) const;

  double fx, fy, cx, cy;
};


class  EdgeStereoSE3ProjectXYZ: public  BaseBinaryEdge<3, Vector3d, VertexSBAPointXYZ, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeStereoSE3ProjectXYZ();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
    const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
    Vector3d obs(_measurement);
    _error = obs - cam_project(v1->estimate().map(v2->estimate()),bf);
  }

  bool isDepthPositive() {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
    const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
    return (v1->estimate().map(v2->estimate()))(2)>0.0;
  }


  virtual void linearizeOplus();

  Vector3d cam_project(const Vector3d & trans_xyz, const double &bf) const;

  double fx, fy, cx, cy, bf;
};

class  EdgeSE3ProjectXYZOnlyPose: public  BaseUnaryEdge<2, Vector2d, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeSE3ProjectXYZOnlyPose(){}

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
    Vector2d obs(_measurement);
    _error = obs-cam_project(v1->estimate().map(Xw));
    //test_value_ += 1;
    //std::cout<<"in compute error test value is "<<test_value_<<",  error[0] is "<<_error(0)<<std::endl;
  }

  bool isDepthPositive() {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
    return (v1->estimate().map(Xw))(2)>0.0;
  }


  virtual void linearizeOplus();

  Vector2d cam_project(const Vector3d & trans_xyz) const;

  Vector3d Xw;
  double fx, fy, cx, cy;

  int test_value_ = 0;
};


class  EdgeStereoSE3ProjectXYZOnlyPose: public  BaseUnaryEdge<3, Vector3d, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeStereoSE3ProjectXYZOnlyPose(){}

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
    Vector3d obs(_measurement);
    _error = obs - cam_project(v1->estimate().map(Xw));
  }

  bool isDepthPositive() {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
    return (v1->estimate().map(Xw))(2)>0.0;
  }


  virtual void linearizeOplus();

  Vector3d cam_project(const Vector3d & trans_xyz) const;

  Vector3d Xw;
  double fx, fy, cx, cy, bf;
};

class  EdgeSE3ProjectIntensityOnlyPoseNID: public  BaseUnaryEdge<1, Eigen::VectorXd, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeSE3ProjectIntensityOnlyPoseNID(){}

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError() {
    ClearPrevH();
    //printf("clear H finish \n");
    ComputeH();

    //if nan value happens, normally it means all the pixel is the last image is mapped to the outside of the current image
    _error(0,0) = (2*H_joint_ - H_ref_ - H_current_)/H_joint_; //-cam_project(v1->estimate().map(Xw));

    // if( isnan(_error(0,0)) )
    //   this->setLevel(1);
    //std::cout<<"H last, current, joint is "<<H_ref_<<","<<H_current_<<","<<H_joint_<<std::endl;
    //std::cout<<"error is "<<_error(0,0)<<std::endl;
    //std::cout<<"mutual information is "<<H_ref_ + H_current_ - H_joint_<<std::endl;
  }

  void set_bspline_relates(std::vector<std::vector<std::vector<float> > >& bs, std::vector<std::vector<std::vector<float> > >& bs_der, int bs_degree, int bin_num);

  void ComputeH();

  void computeHref();

  void ClearPrevH();

  virtual void linearizeOplus();

  double fx_, fy_, cx_, cy_;

  Eigen::VectorXd hg_intensity_current_;
  std::vector<cv::KeyPoint> hg_points_;
  Eigen::Matrix<double, -1, 4> bs_value_current_;
  Eigen::Matrix<double, -1, 4> bs_value_ref_;

  //std::vector<Eigen::Vector3d> x_world_set_;
  std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > x_world_set_;

  cv::Mat image_current_;

  double H_current_ = 0.0, H_ref_ = 0.0, H_joint_ = 0.0;

  float test_H_ref_ = 0.0;

  //b spline coeffcient, b spline derivatives, b spline degree
  std::vector<std::vector<std::vector<float> > > bs_, bs_der_;
  //current intensity probability in each b spline function, last image intensity probability in each b spline function, derivate of last frame to mapped intensity (pixel intensity is mapped from 0~255 to (bin_num - bs_degree) )
  std::vector<float> pro_current_, pro_ref_;
  //joint pobability of b spline function
  std::vector<std::vector<float> > pro_joint_;
  int bs_degree_;
  int bin_num_;

  int ob = 0; //mapped points that are out of image boundary

  double sigma_ = 1e-30;

private:
  // get a gray scale value from reference image (bilinear interpolated)
  inline float get_interpolated_pixel_value ( float x, float y )
  {
    int ix = (int)x;
    int iy = (int)y;
    
    float dx = x - ix;
    float dy = y - iy;
    float dxdy = dx*dy; 

    float xx = x - floor ( x );
    float yy = y - floor ( y );

    return float (
      dxdy * image_current_.ptr<float>(iy+1)[ix+1] 
      + (dy - dxdy) * image_current_.ptr<float>(iy+1)[ix]
      + (dx - dxdy) * image_current_.ptr<float>(iy)[ix+1]
      + (1 - dx - dy + dxdy) * image_current_.ptr<float>(iy)[ix]
    );
  }
};

// class  EdgeSE3ProjectIntensityOnlyPose: public  BaseUnaryEdge<1, Eigen::Matrix<double, 1, 1>, VertexSE3Expmap>{
// public:
//   EIGEN_MAKE_ALIGNED_OPERATOR_NEW

//   EdgeSE3ProjectIntensityOnlyPose(){}

//   bool read(std::istream& is);

//   bool write(std::ostream& os) const;

//   void computeError() {
//     VertexSE3Expmap * vi = static_cast<VertexSE3Expmap *>(_vertices[0]);
//     auto T_lc = vi->estimate().to_homogeneous_matrix();

//     float x_p = hg_depth_current_ * (hg_point_current_.pt.x - cx_) / fx_;
//     float y_p = hg_depth_current_ * (hg_point_current_.pt.y - cy_) / fy_;
//     float z_p = hg_depth_current_;
//     //3d point in pFrame and lFrame
//     Eigen::Matrix<double,4,1> p_p, p_l;
//     p_p << x_p, y_p, z_p, 1;
//     p_l = T_lc * p_p;
//     //std::cout<<T_lc<<std::endl;
//     u_l_ =  p_l(0,0)/p_l(2,0);
//     v_l_ =  p_l(1,0)/p_l(2,0);
//     inv_depth_ = 1.0/p_l(2,0);

//     //2d pixel position in l frame
//     float x = fx_ * p_l(0,0) / p_l(2,0) + cx_;
//     float y = fy_ * p_l(1,0) / p_l(2,0) + cy_;

//     //interpolation, DSO, getInterpolatedElement33() function, bilinear interpolation
//     int ix = (int)x;
//     int iy = (int)y;
    
//     float dx = x - ix;
//     float dy = y - iy;
//     float dxdy = dx*dy; 

//     float intensity_pro = 0.0;
//     //bilinear interporlation of pixel
//     if(ix >= 0 && ix+2 <= image_last_.cols - 1 && iy >= 0 && iy+2 <= image_last_.rows -1){
//       intensity_pro = dxdy * image_last_.ptr<float>(iy+1)[ix+1] 
//                             + (dy - dxdy) * image_last_.ptr<float>(iy+1)[ix]
//                             + (dx - dxdy) * image_last_.ptr<float>(iy)[ix+1]
//                             + (1 - dx - dy + dxdy) * image_last_.ptr<float>(iy)[ix];
//     }
//     else{
//       this->setLevel(1); //else error will be big and not be used
//     }

//     //std::cout<<"original x,y"<<hg_point_current_.pt.x<<", "<<hg_point_current_.pt.y<<", new x and y "<<x<<", "<<y<<std::endl;
//     //std::cout<<"original intensity and current intensity "<<_measurement(0,0)<<", "<<intensity_pro<<std::endl;
    
//     _error(0,0) = intensity_pro - _measurement(0,0);
//     //std::cout<<"error is "<<_error(0,0)<<std::endl;
//   }


//   virtual void linearizeOplus();

//   Vector3d Xw;
//   double fx_, fy_, cx_, cy_;

//   cv::KeyPoint hg_point_current_;
//   float hg_depth_current_;

//   cv::Mat image_last_;

//   double sigma_ = 1e-16;

//   float u_l_, v_l_, inv_depth_;
// };

class EdgeSE3ProjectDirectOnlyPose: public BaseUnaryEdge< 1, double, VertexSE3Expmap>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EdgeSE3ProjectDirectOnlyPose(){}

    virtual void computeError()
    {
        const VertexSE3Expmap* v  =static_cast<const VertexSE3Expmap*> ( _vertices[0] );
        Eigen::Vector3d x_local = v->estimate().map ( x_world_ );
        float x = x_local[0]*fx_/x_local[2] + cx_;
        float y = x_local[1]*fy_/x_local[2] + cy_;
        // check x,y is in the image
        if ( x-4<0 || ( x+4 ) >image_current_.cols || ( y-4 ) <0 || ( y+4 ) >image_current_.rows )
        {
            _error ( 0,0 ) = 0.0;
            this->setLevel ( 1 );
        }
        else
        {
            _error ( 0,0 ) = get_interpolated_pixel_value ( x,y ) - _measurement;
        }
    }

    virtual bool read ( std::istream& in ) {}
    virtual bool write ( std::ostream& out ) const {}

    virtual void linearizeOplus();

    double fx_, fy_, cx_, cy_;

    Vector3d x_world_;

    cv::Mat image_current_;

private:
    // get a gray scale value from reference image (bilinear interpolated)
    inline float get_interpolated_pixel_value ( float x, float y )
    {
      int ix = (int)x;
      int iy = (int)y;
      
      float dx = x - ix;
      float dy = y - iy;
      float dxdy = dx*dy; 

      float xx = x - floor ( x );
      float yy = y - floor ( y );

      return float (
        dxdy * image_current_.ptr<float>(iy+1)[ix+1] 
        + (dy - dxdy) * image_current_.ptr<float>(iy+1)[ix]
        + (dx - dxdy) * image_current_.ptr<float>(iy)[ix+1]
        + (1 - dx - dy + dxdy) * image_current_.ptr<float>(iy)[ix]
      );
    }
};


} // end namespace

#endif
