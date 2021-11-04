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

#include "types_six_dof_expmap.h"

#include "../core/factory.h"
#include "../stuff/macros.h"

namespace g2o {

using namespace std;


Vector2d project2d(const Vector3d& v)  {
  Vector2d res;
  res(0) = v(0)/v(2);
  res(1) = v(1)/v(2);
  return res;
}

Vector3d unproject2d(const Vector2d& v)  {
  Vector3d res;
  res(0) = v(0);
  res(1) = v(1);
  res(2) = 1;
  return res;
}

VertexSE3Expmap::VertexSE3Expmap() : BaseVertex<6, SE3Quat>() {
}

bool VertexSE3Expmap::read(std::istream& is) {
  Vector7d est;
  for (int i=0; i<7; i++)
    is  >> est[i];
  SE3Quat cam2world;
  cam2world.fromVector(est);
  setEstimate(cam2world.inverse());
  return true;
}

bool VertexSE3Expmap::write(std::ostream& os) const {
  SE3Quat cam2world(estimate().inverse());
  for (int i=0; i<7; i++)
    os << cam2world[i] << " ";
  return os.good();
}


EdgeSE3ProjectXYZ::EdgeSE3ProjectXYZ() : BaseBinaryEdge<2, Vector2d, VertexSBAPointXYZ, VertexSE3Expmap>() {
}

bool EdgeSE3ProjectXYZ::read(std::istream& is){
  for (int i=0; i<2; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeSE3ProjectXYZ::write(std::ostream& os) const {

  for (int i=0; i<2; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}


void EdgeSE3ProjectXYZ::linearizeOplus() {
  VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[1]);
  SE3Quat T(vj->estimate());
  VertexSBAPointXYZ* vi = static_cast<VertexSBAPointXYZ*>(_vertices[0]);
  Vector3d xyz = vi->estimate();
  Vector3d xyz_trans = T.map(xyz);

  double x = xyz_trans[0];
  double y = xyz_trans[1];
  double z = xyz_trans[2];
  double z_2 = z*z;

  Matrix<double,2,3> tmp;
  tmp(0,0) = fx;
  tmp(0,1) = 0;
  tmp(0,2) = -x/z*fx;

  tmp(1,0) = 0;
  tmp(1,1) = fy;
  tmp(1,2) = -y/z*fy;

  _jacobianOplusXi =  -1./z * tmp * T.rotation().toRotationMatrix();

  _jacobianOplusXj(0,0) =  x*y/z_2 *fx;
  _jacobianOplusXj(0,1) = -(1+(x*x/z_2)) *fx;
  _jacobianOplusXj(0,2) = y/z *fx;
  _jacobianOplusXj(0,3) = -1./z *fx;
  _jacobianOplusXj(0,4) = 0;
  _jacobianOplusXj(0,5) = x/z_2 *fx;

  _jacobianOplusXj(1,0) = (1+y*y/z_2) *fy;
  _jacobianOplusXj(1,1) = -x*y/z_2 *fy;
  _jacobianOplusXj(1,2) = -x/z *fy;
  _jacobianOplusXj(1,3) = 0;
  _jacobianOplusXj(1,4) = -1./z *fy;
  _jacobianOplusXj(1,5) = y/z_2 *fy;
}

Vector2d EdgeSE3ProjectXYZ::cam_project(const Vector3d & trans_xyz) const{
  Vector2d proj = project2d(trans_xyz);
  Vector2d res;
  res[0] = proj[0]*fx + cx;
  res[1] = proj[1]*fy + cy;
  return res;
}


Vector3d EdgeStereoSE3ProjectXYZ::cam_project(const Vector3d & trans_xyz, const double &bf) const{
  const double invz = 1.0f/trans_xyz[2];
  Vector3d res;
  res[0] = trans_xyz[0]*invz*fx + cx;
  res[1] = trans_xyz[1]*invz*fy + cy;
  res[2] = res[0] - bf*invz;
  return res;
}

EdgeStereoSE3ProjectXYZ::EdgeStereoSE3ProjectXYZ() : BaseBinaryEdge<3, Vector3d, VertexSBAPointXYZ, VertexSE3Expmap>() {
}

bool EdgeStereoSE3ProjectXYZ::read(std::istream& is){
  for (int i=0; i<=3; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeStereoSE3ProjectXYZ::write(std::ostream& os) const {

  for (int i=0; i<=3; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}

void EdgeStereoSE3ProjectXYZ::linearizeOplus() {
  VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[1]);
  SE3Quat T(vj->estimate());
  VertexSBAPointXYZ* vi = static_cast<VertexSBAPointXYZ*>(_vertices[0]);
  Vector3d xyz = vi->estimate();
  Vector3d xyz_trans = T.map(xyz);

  const Matrix3d R =  T.rotation().toRotationMatrix();

  double x = xyz_trans[0];
  double y = xyz_trans[1];
  double z = xyz_trans[2];
  double z_2 = z*z;

  _jacobianOplusXi(0,0) = -fx*R(0,0)/z+fx*x*R(2,0)/z_2;
  _jacobianOplusXi(0,1) = -fx*R(0,1)/z+fx*x*R(2,1)/z_2;
  _jacobianOplusXi(0,2) = -fx*R(0,2)/z+fx*x*R(2,2)/z_2;

  _jacobianOplusXi(1,0) = -fy*R(1,0)/z+fy*y*R(2,0)/z_2;
  _jacobianOplusXi(1,1) = -fy*R(1,1)/z+fy*y*R(2,1)/z_2;
  _jacobianOplusXi(1,2) = -fy*R(1,2)/z+fy*y*R(2,2)/z_2;

  _jacobianOplusXi(2,0) = _jacobianOplusXi(0,0)-bf*R(2,0)/z_2;
  _jacobianOplusXi(2,1) = _jacobianOplusXi(0,1)-bf*R(2,1)/z_2;
  _jacobianOplusXi(2,2) = _jacobianOplusXi(0,2)-bf*R(2,2)/z_2;

  _jacobianOplusXj(0,0) =  x*y/z_2 *fx;
  _jacobianOplusXj(0,1) = -(1+(x*x/z_2)) *fx;
  _jacobianOplusXj(0,2) = y/z *fx;
  _jacobianOplusXj(0,3) = -1./z *fx;
  _jacobianOplusXj(0,4) = 0;
  _jacobianOplusXj(0,5) = x/z_2 *fx;

  _jacobianOplusXj(1,0) = (1+y*y/z_2) *fy;
  _jacobianOplusXj(1,1) = -x*y/z_2 *fy;
  _jacobianOplusXj(1,2) = -x/z *fy;
  _jacobianOplusXj(1,3) = 0;
  _jacobianOplusXj(1,4) = -1./z *fy;
  _jacobianOplusXj(1,5) = y/z_2 *fy;

  _jacobianOplusXj(2,0) = _jacobianOplusXj(0,0)-bf*y/z_2;
  _jacobianOplusXj(2,1) = _jacobianOplusXj(0,1)+bf*x/z_2;
  _jacobianOplusXj(2,2) = _jacobianOplusXj(0,2);
  _jacobianOplusXj(2,3) = _jacobianOplusXj(0,3);
  _jacobianOplusXj(2,4) = 0;
  _jacobianOplusXj(2,5) = _jacobianOplusXj(0,5)-bf/z_2;
}


//Only Pose

bool EdgeSE3ProjectXYZOnlyPose::read(std::istream& is){
  for (int i=0; i<2; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeSE3ProjectXYZOnlyPose::write(std::ostream& os) const {

  for (int i=0; i<2; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}


void EdgeSE3ProjectXYZOnlyPose::linearizeOplus() {
  VertexSE3Expmap * vi = static_cast<VertexSE3Expmap *>(_vertices[0]);
  Vector3d xyz_trans = vi->estimate().map(Xw);

  //std::cout<<"to be estimated matrix \n"<<vi->estimate().to_homogeneous_matrix()<<std::endl;

  double x = xyz_trans[0];
  double y = xyz_trans[1];
  double invz = 1.0/xyz_trans[2];
  double invz_2 = invz*invz;

  _jacobianOplusXi(0,0) =  x*y*invz_2 *fx;
  _jacobianOplusXi(0,1) = -(1+(x*x*invz_2)) *fx;
  _jacobianOplusXi(0,2) = y*invz *fx;
  _jacobianOplusXi(0,3) = -invz *fx;
  _jacobianOplusXi(0,4) = 0;
  _jacobianOplusXi(0,5) = x*invz_2 *fx;

  _jacobianOplusXi(1,0) = (1+y*y*invz_2) *fy;
  _jacobianOplusXi(1,1) = -x*y*invz_2 *fy;
  _jacobianOplusXi(1,2) = -x*invz *fy;
  _jacobianOplusXi(1,3) = 0;
  _jacobianOplusXi(1,4) = -invz *fy;
  _jacobianOplusXi(1,5) = y*invz_2 *fy;

  //std::cout<<"the jacobian value, x,y,z only pose is \n"<<_jacobianOplusXi<<std::endl;
}

Vector2d EdgeSE3ProjectXYZOnlyPose::cam_project(const Vector3d & trans_xyz) const{
  Vector2d proj = project2d(trans_xyz);
  Vector2d res;
  res[0] = proj[0]*fx + cx;
  res[1] = proj[1]*fy + cy;
  return res;
}


Vector3d EdgeStereoSE3ProjectXYZOnlyPose::cam_project(const Vector3d & trans_xyz) const{
  const double invz = 1.0f/trans_xyz[2];
  Vector3d res;
  res[0] = trans_xyz[0]*invz*fx + cx;
  res[1] = trans_xyz[1]*invz*fy + cy;
  res[2] = res[0] - bf*invz;
  return res;
}


bool EdgeStereoSE3ProjectXYZOnlyPose::read(std::istream& is){
  for (int i=0; i<=3; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeStereoSE3ProjectXYZOnlyPose::write(std::ostream& os) const {

  for (int i=0; i<=3; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}

void EdgeStereoSE3ProjectXYZOnlyPose::linearizeOplus() {
  VertexSE3Expmap * vi = static_cast<VertexSE3Expmap *>(_vertices[0]);
  Vector3d xyz_trans = vi->estimate().map(Xw);

  double x = xyz_trans[0];
  double y = xyz_trans[1];
  double invz = 1.0/xyz_trans[2];
  double invz_2 = invz*invz;

  _jacobianOplusXi(0,0) =  x*y*invz_2 *fx;
  _jacobianOplusXi(0,1) = -(1+(x*x*invz_2)) *fx;
  _jacobianOplusXi(0,2) = y*invz *fx;
  _jacobianOplusXi(0,3) = -invz *fx;
  _jacobianOplusXi(0,4) = 0;
  _jacobianOplusXi(0,5) = x*invz_2 *fx;

  _jacobianOplusXi(1,0) = (1+y*y*invz_2) *fy;
  _jacobianOplusXi(1,1) = -x*y*invz_2 *fy;
  _jacobianOplusXi(1,2) = -x*invz *fy;
  _jacobianOplusXi(1,3) = 0;
  _jacobianOplusXi(1,4) = -invz *fy;
  _jacobianOplusXi(1,5) = y*invz_2 *fy;

  _jacobianOplusXi(2,0) = _jacobianOplusXi(0,0)-bf*y*invz_2;
  _jacobianOplusXi(2,1) = _jacobianOplusXi(0,1)+bf*x*invz_2;
  _jacobianOplusXi(2,2) = _jacobianOplusXi(0,2);
  _jacobianOplusXi(2,3) = _jacobianOplusXi(0,3);
  _jacobianOplusXi(2,4) = 0;
  _jacobianOplusXi(2,5) = _jacobianOplusXi(0,5)-bf*invz_2;
}

/***Add by zhaozhong chen**/
bool EdgeSE3ProjectIntensityOnlyPoseNID::read(std::istream& is){
  // for (int i=0; i<2; i++){
  //   is >> _measurement[i];
  // }
  // for (int i=0; i<2; i++)
  //   for (int j=i; j<2; j++) {
  //     is >> information()(i,j);
  //     if (i!=j)
  //       information()(j,i)=information()(i,j);
  //   }
  return true;
}

bool EdgeSE3ProjectIntensityOnlyPoseNID::write(std::ostream& os) const {

  // for (int i=0; i<2; i++){
  //   os << measurement()[i] << " ";
  // }

  // for (int i=0; i<2; i++)
  //   for (int j=i; j<2; j++){
  //     os << " " <<  information()(i,j);
  //   }
  return os.good();
}


void EdgeSE3ProjectIntensityOnlyPoseNID::linearizeOplus() {
  //std::cout<<"in g2o compute jacobian.............."<<std::endl;
  VertexSE3Expmap * vi = static_cast<VertexSE3Expmap *>(_vertices[0]);
  Eigen::VectorXd obs(_measurement);
  
  std::vector<std::vector<float> > d_sum_bs_pose(bin_num_, std::vector<float>(6, 0.0));

  //bin_num * bin_num * 6 dimension
  std::vector<std::vector<std::vector<float> > > d_sum_joint_bs_pose(bin_num_, std::vector<std::vector<float> >(bin_num_, std::vector<float>(6,0.0)));

  //derivative, mappde intensity to intensity to 
  float d_mi_i = (bin_num_- bs_degree_)/255.0;

  for(int i = 0 ; i<_measurement.rows(); i++){

    Eigen::Vector3d p_c = vi->estimate().map ( x_world_set_[i] );

    if(obs(i,0) >= 255)
      obs(i,0) = 254.999;
    if(obs(i,0) < 0)
      obs(i, 0) = 0.0;
    
    float bin_pos_ref =  obs(i,0)* (bin_num_- bs_degree_)/255.0;
    int bins_index_ref = floor(bin_pos_ref);

    float u_c =  p_c(0,0)/p_c(2,0);
    float v_c =  p_c(1,0)/p_c(2,0);

    // if(p_l(2.0) == 0){
    //   std::cout<<"warning "<<std::endl;
    // }

    double x = p_c(0,0);
    double y = p_c(1,0);
    double invz = 1.0/p_c(2,0);
    double invz_2 = invz*invz;

    // jacobian from se3 to u,v
    // NOTE that in g2o the Lie algebra is (\omega, \epsilon), where \omega is so(3) and \epsilon the translation
    Eigen::Matrix<double, 2, 6> jacobian_uv_ksai;


    //2d pixel position in l frame
    float u = fx_ * u_c + cx_;
    float v = fy_ * v_c + cy_;

    float bin_pos_current =  hg_intensity_current_(i,0)* (bin_num_-3.0)/255.0;
    float pos_cubic_current = bin_pos_current * bin_pos_current * bin_pos_current;
    float pos_qua_current  = bin_pos_current * bin_pos_current;
    int bins_index_current = floor(bin_pos_current);

    float gradient_x, gradient_y;

    Eigen::Matrix<double, 1, 2> jacobian_pixel_uv;
    //calculate the four corner pixel's gradient and do bilinear inerpolation
    if(u >= 0 && u+3 <= image_current_.cols - 1 && v >= 0 && v+3 <= image_current_.rows){
      jacobian_pixel_uv ( 0,0 ) = ( get_interpolated_pixel_value ( u+1,v )-get_interpolated_pixel_value ( u-1,v ) ) /2;
      jacobian_pixel_uv ( 0,1 ) = ( get_interpolated_pixel_value ( u,v+1 )-get_interpolated_pixel_value ( u,v-1 ) ) /2;
      //std::cout<<"gradient x and y "<<gradient_x<<", "<<gradient_y<<std::endl;

      jacobian_uv_ksai ( 0,0 ) = - x*y*invz_2 *fx_;
      jacobian_uv_ksai ( 0,1 ) = ( 1+ ( x*x*invz_2 ) ) *fx_;
      jacobian_uv_ksai ( 0,2 ) = - y*invz *fx_;
      jacobian_uv_ksai ( 0,3 ) = invz *fx_;
      jacobian_uv_ksai ( 0,4 ) = 0;
      jacobian_uv_ksai ( 0,5 ) = -x*invz_2 *fx_;

      jacobian_uv_ksai ( 1,0 ) = - ( 1+y*y*invz_2 ) *fy_;
      jacobian_uv_ksai ( 1,1 ) = x*y*invz_2 *fy_;
      jacobian_uv_ksai ( 1,2 ) = x*invz *fy_;
      jacobian_uv_ksai ( 1,3 ) = 0;
      jacobian_uv_ksai ( 1,4 ) = invz *fy_;
      jacobian_uv_ksai ( 1,5 ) = -y*invz_2 *fy_;
    }
    else{
      jacobian_pixel_uv ( 0,0 ) = 0.0;
      jacobian_pixel_uv ( 0,1 ) = 0.0;
      jacobian_uv_ksai = Eigen::Matrix<double, 2, 6>::Zero();
      //may need think more about this
      //std::cout<<"opps..... out of image boundary"<<std::endl;
      continue;
    }
    
    Eigen::Matrix<double, 1, 6> d_i_pose;
    d_i_pose = jacobian_pixel_uv * jacobian_uv_ksai;

    //std::cout<<"intensity to pose jacobian is "<<d_i_pose<<std::endl;

    Eigen::MatrixXd sum_d_i_pose = Eigen::MatrixXd::Zero(6,1);//sum of derivatives of intensity to pose

    // std::cout<<"gradient pixel to pose ";
    // for(int i = 0; i<6; i++){
    //   std::cout<<d_i_pose[i]<<",";
    // }
    // std::cout<<std::endl;

    std::vector<std::vector<float> > d_coe0_current, d_coe1_current, d_coe2_current, d_coe3_current;


    d_coe0_current = bs_der_[bins_index_current];
    d_coe1_current = bs_der_[bins_index_current + 1];
    d_coe2_current = bs_der_[bins_index_current + 2];
    d_coe3_current = bs_der_[bins_index_current + 3];

    std::vector<float> d_bs_mi(bin_num_, 0.0);

    //Checked when bin position in (0,1) (2,3) (3,4), all correct
    //std::cout<<"measurement is "<<hg_intensity_last_(i,0)<<", bin position is "<<bin_pos_last<<std::endl;
    d_bs_mi[bins_index_current]     = d_coe0_current[3][0] * pos_qua_current + d_coe0_current[3][1] * bin_pos_current + d_coe0_current[3][2];
    d_bs_mi[bins_index_current + 1] = d_coe1_current[2][0] * pos_qua_current + d_coe1_current[2][1] * bin_pos_current + d_coe1_current[2][2];
    d_bs_mi[bins_index_current + 2] = d_coe2_current[1][0] * pos_qua_current + d_coe2_current[1][1] * bin_pos_current + d_coe2_current[1][2];
    d_bs_mi[bins_index_current + 3] = d_coe3_current[0][0] * pos_qua_current + d_coe3_current[0][1] * bin_pos_current + d_coe3_current[0][2];
    // std::cout<<"bs derivative value "<<d_bs_mi[bins_index_last]<<", "<<d_bs_mi[bins_index_last+1]<<", "<<d_bs_mi[bins_index_last+2]<<", "<<d_bs_mi[bins_index_last+3]<<std::endl;

    for(int m = 0; m<bs_degree_+1; m++)
      for(int n = 0; n<6; n++){
        d_sum_bs_pose[bins_index_current + m][n] += d_bs_mi[bins_index_current + m] * d_mi_i * d_i_pose(0,n);
      }

    for(int k = 0; k<bs_degree_+1; k++)
      for(int m = 0; m<bs_degree_+1; m++)
        for(int n = 0; n<6; n++){
          d_sum_joint_bs_pose[bins_index_ref + k][bins_index_current + m][n] += bs_value_ref_(i,k) * d_bs_mi[bins_index_current + m] * d_mi_i * d_i_pose[n];
        }
  }

  //pro_last_.size() = bin_num
  // for(int i = 0; i < bin_num_ ; i++){
  //   pro_last_[i] /= _measurement.rows();
  // }

  // for(int i = 0; i < bin_num_; i++)
  //   for(int j = 0; j<bin_num_; j++){
  //     pro_joint_[i][j] /= _measurement.rows(); 
  //   }
  
  // for(int i = 0; i < bin_num_ ; i++){
  //   H_last_ -= pro_last_[i] * log2(pro_last_[i]);
  // }

  // for(int i = 0; i < bin_num_; i++)
  //   for(int j = 0; j<bin_num_; j++){
  //     H_joint_ -= pro_joint_[i][j] * log2(pro_joint_[i][j]);
  //   }

  for(int m = 0; m<bin_num_; m++)
    for(int n = 0; n<6; n++){
      d_sum_bs_pose[m][n] /= (_measurement.rows() - ob);
      //std::cout<<"sum bspline to pose derivarive"<<d_sum_bs_pose[m][n]<<" for bin num "<<m<<" derivative num "<<n<<std::endl;
  }

  for(int m = 0; m<bin_num_; m++)
    for(int n = 0; n<bin_num_; n++)
      for(int k = 0; k<6; k++){
        d_sum_joint_bs_pose[m][n][k] /= (_measurement.rows() - ob);
        //std::cout<<"joint sum bspline to pose derivarive"<<d_sum_joint_bs_pose[m][n][k]<<" for bin num "<<m<<","<<n<<", der number "<<k<<std::endl;
      }

  //derivative H_joint to pose
  std::vector<double> d_hj_p(6, 0.0);
  for(int i = 0; i < 6; i++){
    float tmp = 0.0;
    for(int m = 0; m<bin_num_; m++){
      for(int n = 0; n<bin_num_; n++){
        if(pro_joint_[m][n] < sigma_)
          continue;
        tmp -= (1.0 + log2(pro_joint_[m][n])) * d_sum_joint_bs_pose[m][n][i];
      }
    }
    d_hj_p[i] = tmp;
  }  

  //derivative H_last_ to pose
  std::vector<double> d_hl_p(6, 0.0);
  for(int i = 0; i < 6; i++){
    for(int j = 0; j < bin_num_; j++){
      if(pro_current_[j] < sigma_)
        continue;
      d_hl_p[i] -= (1.0 + log2(pro_current_[j])) * d_sum_bs_pose[j][i];
    }
  }

  //std::cout<<"joint derivative "<<d_hj_p[0]<<", "<<d_hj_p[1]<<", "<<d_hj_p[2]<<d_hj_p[3]<<", "<<d_hj_p[4]<<", "<<d_hj_p[5]<<std::endl;
  //std::cout<<"last derivative "<<d_hl_p[0]<<", "<<d_hl_p[1]<<", "<<d_hl_p[2]<<d_hl_p[3]<<", "<<d_hl_p[4]<<", "<<d_hl_p[5]<<std::endl;
  //std::cout<<"H last current joint  "<<H_ref_<<", "<<H_current_<<","<<H_joint_<<std::endl;
  float inv_square_hj = 1.0/(H_joint_ * H_joint_);
  _jacobianOplusXi(0,0) = (d_hj_p[0] * (H_current_ + H_ref_) - d_hl_p[0] * H_joint_) * inv_square_hj;
  _jacobianOplusXi(0,1) = (d_hj_p[1] * (H_current_ + H_ref_) - d_hl_p[1] * H_joint_) * inv_square_hj;
  _jacobianOplusXi(0,2) = (d_hj_p[2] * (H_current_ + H_ref_) - d_hl_p[2] * H_joint_) * inv_square_hj; 
  _jacobianOplusXi(0,3) = (d_hj_p[3] * (H_current_ + H_ref_) - d_hl_p[3] * H_joint_) * inv_square_hj; 
  _jacobianOplusXi(0,4) = (d_hj_p[4] * (H_current_ + H_ref_) - d_hl_p[4] * H_joint_) * inv_square_hj; 
  _jacobianOplusXi(0,5) = (d_hj_p[5] * (H_current_ + H_ref_) - d_hl_p[5] * H_joint_) * inv_square_hj;

  //_jacobianOplusXi = 1000 * _jacobianOplusXi;   

  //std::cout<<"the jacobian is \n"<<_jacobianOplusXi<<std::endl;
  // if(isnan(_jacobianOplusXi(0,0)) || isnan(_jacobianOplusXi(0,1)) || isnan(_jacobianOplusXi(0,2)) || isnan(_jacobianOplusXi(0,3)) || isnan(_jacobianOplusXi(0,4)) || isnan(_jacobianOplusXi(0,5))){
  //   std::cout<<"exit in jacobian estimation "<<std::endl;
  //   std::cout<<"the matrix to be estimated is \n"<<vi->estimate().to_homogeneous_matrix()<<std::endl;
  //   std::cout<<"the Href, Hcurrent, Hjoint is "<<H_ref_<<","<<H_current_<<","<<H_joint_<<std::endl;
  //   exit(0);
  // } 
  

}


void EdgeSE3ProjectIntensityOnlyPoseNID::ComputeH(){

  const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
  Eigen::VectorXd obs(_measurement);
  ob = 0;
  //std::cout<<"mapped maptrix is \n"<<v1->estimate().to_homogeneous_matrix()<<std::endl;

  std::vector<float> test_pro_current(bin_num_,0);
  std::vector<std::vector<float> > test_pro_joint(bin_num_, std::vector<float>(bin_num_,0));
  float test_H_current = 0.0, test_H_joint = 0.0;
  for(int i = 0 ; i<obs.rows(); i++){
    //current frame mutual information probabilit
    Eigen::Vector3d p_c = v1->estimate().map ( x_world_set_[i] );

    if(obs(i,0) >= 255)
      obs(i,0) = 254.999;
    if(obs(i,0) < 0)
      obs(i, 0) = 0.0;

    float bin_pos_ref =  obs(i,0)* (bin_num_- bs_degree_)/255.0;
    int bins_index_ref = floor(bin_pos_ref);
    
    //2d pixel position in current frame
    float u = fx_ * p_c(0,0) / p_c(2,0) + cx_;
    float v = fy_ * p_c(1,0) / p_c(2,0) + cy_;

    //if(i%5 == 0){
      //std::cout<<"original 3d point is "<<x_world_set_[i].transpose()<<std::endl;
      //std::cout<<"the mapped point 3d is "<<p_c.transpose()<<std::endl;
      //std::cout<<"the corresponding 2d point in current frame is "<<hg_points_[i].pt.x<<", "<<hg_points_[i].pt.y<<std::endl;
      //std::cout<<"aftre mapping, the pixel x and y is "<<u<<","<<v<<std::endl;
    //}

    //bilinear interporlation of pixel. DSO getInterpolatedElement33() function, bilinear interpolation
    if(u >= 0 && u+3 <= image_current_.cols && v >= 0 && v+3 <= image_current_.rows){
      hg_intensity_current_(i,0) = get_interpolated_pixel_value(u,v);
    }
    else{
      ob++;
      continue;
    }

    // if(iy+1>image_last_.rows-1 || iy < 0 || ix>image_last_.cols-1 || ix < 0){
    //   std::cout<<" the intensity is  "<< hg_intensity_last_(i,0)<<std::endl;
    // }
    if(hg_intensity_current_(i,0) >= 255)
      hg_intensity_current_(i,0) = 254.999;
    if(hg_intensity_current_(i,0) < 0)
      hg_intensity_current_(i, 0) = 0.0;

    //std::cout<<"original point, mapped point, original intensity, mapped intensity"
            // <<hg_points_[i].pt.x<<", "<<hg_points_[i].pt.y<<",,,"<<u<<", "<<v<<",,,"<<obs(i,0)<<","<<hg_intensity_current_(i,0)<<std::endl;
    //last frame mutual information probability
    float bin_pos_current = hg_intensity_current_(i,0)* (bin_num_-3.0)/255.0;
    float pos_cubic_current = bin_pos_current * bin_pos_current * bin_pos_current;
    float pos_qua_current  = bin_pos_current * bin_pos_current;
    float bins_index_current = floor(bin_pos_current);
    
    std::vector<std::vector<float> > coe0_current, coe1_current, coe2_current, coe3_current;//coefficient of each basis function and its derivatives

    coe0_current = bs_[bins_index_current];
    coe1_current = bs_[bins_index_current + 1];
    coe2_current = bs_[bins_index_current + 2];
    coe3_current = bs_[bins_index_current + 3];

    float bs_value0_current = coe0_current[3][0] * pos_cubic_current + coe0_current[3][1] * pos_qua_current + coe0_current[3][2] * bin_pos_current + coe0_current[3][3];
    bs_value_current_(i, 0 ) = bs_value0_current;
    float bs_value1_current = coe1_current[2][0] * pos_cubic_current + coe1_current[2][1] * pos_qua_current + coe1_current[2][2] * bin_pos_current + coe1_current[2][3];
    bs_value_current_(i, 1 ) = bs_value1_current;
    float bs_value2_current = coe2_current[1][0] * pos_cubic_current + coe2_current[1][1] * pos_qua_current + coe2_current[1][2] * bin_pos_current + coe2_current[1][3];
    bs_value_current_(i, 2) = bs_value2_current;
    float bs_value3_current = coe3_current[0][0] * pos_cubic_current + coe3_current[0][1] * pos_qua_current + coe3_current[0][2] * bin_pos_current + coe3_current[0][3];
    bs_value_current_(i, 3) = bs_value3_current;

    pro_current_[bins_index_current]   += bs_value0_current;
    pro_current_[bins_index_current+1] += bs_value1_current;
    pro_current_[bins_index_current+2] += bs_value2_current;
    pro_current_[bins_index_current+3] += bs_value3_current;
    test_pro_current[bins_index_current]  += 1.0;

    // std::cout<<"bs value "<<bs_value0_current<<", "<<bs_value1_current<<","<<bs_value2_current<<","<<bs_value3_current<<std::endl;
    // std::cout<<"pro current "<<pro_current_[bins_index_current]<<","<<pro_current_[bins_index_current+1]<<","<<pro_current_[bins_index_current+2]<<","<<pro_current_[bins_index_current+3]<<std::endl;
    
    //(bs_degree + 1) * (bs_degree + 1) combination
    for(int m = 0 ; m<bs_degree_+1; m++)
      for(int n = 0; n<bs_degree_+1; n++){
        pro_joint_[bins_index_ref+m][bins_index_current+n] += bs_value_ref_(i,m)*bs_value_current_(i,n);
      }
    test_pro_joint[bins_index_ref][bins_index_current] += 1.0;
    
  }

  //std::cout<<"final ob is "<<ob<<std::endl;

  //if too many points are out of image boundary, we don't use this edge
  if(obs.rows() - ob < 300){
    //std::cout<<"no enough points in the image boundary"<<std::endl;
    this->setLevel(1);
  }

  //std::cout<<"ob is  "<<ob<<", row number is "<<obs.rows()<<std::endl;
  // std::cout<<"value of procurrent after";
  // for(int i = 0; i < bin_num_ ; i++){
  //   std::cout<<pro_current_[i]<<",";
  // }
  // std::cout<<std::endl;
  

  //pro_last_.size() = bin_num
  for(int i = 0; i < bin_num_ ; i++){
    pro_current_[i] /= (obs.rows() - ob);
  }

  for(int i = 0; i < bin_num_; i++)
    for(int j = 0; j<bin_num_; j++){
      pro_joint_[i][j] /= (obs.rows() - ob); 
    }
  
  for(int i = 0; i < bin_num_ ; i++){
    if(pro_current_[i] < sigma_)
      continue;
    H_current_ -= pro_current_[i] * log2(pro_current_[i]);
  }

  for(int i = 0; i < bin_num_; i++)
    for(int j = 0; j<bin_num_; j++){
      if(pro_joint_[i][j] < sigma_)
        continue;
      H_joint_ -= pro_joint_[i][j] * log2(pro_joint_[i][j]);
    }

  //TEST
  for(int i = 0; i < bin_num_ ; i++){
    test_pro_current[i] /= (obs.rows() - ob);
  }

  for(int i = 0; i < bin_num_; i++)
    for(int j = 0; j<bin_num_; j++){
      test_pro_joint[i][j] /= (obs.rows() - ob); 
    }
  
  for(int i = 0; i < bin_num_ ; i++){
    if(test_pro_current[i] < sigma_)
      continue;
    test_H_current -= test_pro_current[i] * log2(test_pro_current[i]);
  }

  for(int i = 0; i < bin_num_; i++)
    for(int j = 0; j<bin_num_; j++){
      if(test_pro_joint[i][j] < sigma_)
        continue;
      test_H_joint -= test_pro_joint[i][j] * log2(test_pro_joint[i][j]);
    }

  // std::cout<<"Href, current, joint from B spline "<<H_ref_<<","<<H_current_<<","<<H_joint_<<", MI "<<H_ref_ + H_current_ - H_joint_<<", NID "<<(2*H_joint_ - H_ref_ - H_current_)/H_joint_<<std::endl;
  // std::cout<<"Href, current, joint from standard approach "<<test_H_ref_<<","<<test_H_current<<","<<test_H_joint<<", MI "<<test_H_current+test_H_ref_-test_H_joint<<", NID "<< (2*test_H_joint - test_H_current - test_H_ref_)/test_H_joint<<std::endl;
  //std::cout<<"\n"<<std::endl;
}

void EdgeSE3ProjectIntensityOnlyPoseNID::set_bspline_relates(std::vector<std::vector<std::vector<float> > >& bs, std::vector<std::vector<std::vector<float> > >& bs_der, int bs_degree, int bin_num){
  bs_ = bs;
  bs_der_ = bs_der;
  bs_degree_ = bs_degree;
  bin_num_ = bin_num;

  pro_current_ = std::vector<float>(bin_num_,0.0);
  pro_ref_ = std::vector<float>(bin_num_,0.0);
  pro_joint_ = std::vector<std::vector<float> >(bin_num_);
  for(int i = 0; i < bin_num_; i++){
    pro_joint_[i] = std::vector<float>(bin_num_, 0.0);
  }

  bs_value_ref_ = Eigen::MatrixXd::Zero(_measurement.rows(),4);
  bs_value_current_ = Eigen::MatrixXd::Zero(_measurement.rows(),4);
  hg_intensity_current_ = Eigen::VectorXd::Zero(_measurement.rows());
};

void EdgeSE3ProjectIntensityOnlyPoseNID::computeHref(){
  Eigen::VectorXd obs(_measurement);

  std::vector<float> test_pro_ref(bin_num_,0);
  
  for(int i = 0 ; i<obs.rows(); i++){
    //current frame mutual information probability
    if(obs(i,0) >= 255)
      obs(i,0) = 254.999;
    if(obs(i,0) < 0)
      obs(i, 0) = 0.0;
    float bin_pos_ref =  obs(i,0) * (bin_num_- bs_degree_)/255.0;
    float pos_cubic_ref = bin_pos_ref * bin_pos_ref * bin_pos_ref;
    float pos_qua_ref  = bin_pos_ref * bin_pos_ref;
    int bins_index_ref = floor(bin_pos_ref);
    std::vector<std::vector<float> > coe0_ref, coe1_ref, coe2_ref, coe3_ref;//coefficient of each basis function and its derivatives
    coe0_ref = bs_[bins_index_ref];
    coe1_ref = bs_[bins_index_ref + 1];
    coe2_ref = bs_[bins_index_ref + 2];
    coe3_ref = bs_[bins_index_ref + 3];

    float bs_value0_ref = coe0_ref[3][0] * pos_cubic_ref + coe0_ref[3][1] * pos_qua_ref + coe0_ref[3][2] * bin_pos_ref + coe0_ref[3][3];
    bs_value_ref_(i,0) = bs_value0_ref;
    float bs_value1_ref = coe1_ref[2][0] * pos_cubic_ref + coe1_ref[2][1] * pos_qua_ref + coe1_ref[2][2] * bin_pos_ref + coe1_ref[2][3];
    bs_value_ref_(i,1) = bs_value1_ref;
    float bs_value2_ref = coe2_ref[1][0] * pos_cubic_ref + coe2_ref[1][1] * pos_qua_ref + coe2_ref[1][2] * bin_pos_ref + coe2_ref[1][3];
    bs_value_ref_(i,2) = bs_value2_ref;
    float bs_value3_ref = coe3_ref[0][0] * pos_cubic_ref + coe3_ref[0][1] * pos_qua_ref + coe3_ref[0][2] * bin_pos_ref + coe3_ref[0][3];
    bs_value_ref_(i,3) = bs_value3_ref;

    pro_ref_[bins_index_ref]   += bs_value0_ref;
    pro_ref_[bins_index_ref+1] += bs_value1_ref;
    pro_ref_[bins_index_ref+2] += bs_value2_ref;
    pro_ref_[bins_index_ref+3] += bs_value3_ref;

    test_pro_ref[bins_index_ref] += 1.0;

    // have checked the bs value interpolation matches the number from https://wrf.ecse.rpi.edu/wiki/ComputerGraphicsFall2013/guha/Code/bSplines.cpp 
    //std::cout<<"observation is "<<obs(i,0)<<", bin position is "<<bin_pos_ref<<std::endl;
    //std::cout<<"bs value "<<bs_value0_ref<<", "<<bs_value1_ref<<","<<bs_value2_ref<<", "<<bs_value3_ref<<std::endl;
    //std::cout<<"i "<<i<<std::endl;
  }

  for(int i = 0; i < bin_num_ ; i++)
    pro_ref_[i] /=  obs.rows();
  
  for(int i = 0; i < bin_num_ ; i++){
    if(pro_ref_[i] < sigma_)
      continue;
    H_ref_ -= pro_ref_[i] * log2(pro_ref_[i]);
  }

  for(int i = 0; i < bin_num_ ; i++)
    test_pro_ref[i] /=  obs.rows();
  
  for(int i = 0; i < bin_num_ ; i++){
    if(test_pro_ref[i] < sigma_)
      continue;
    test_H_ref_ -= test_pro_ref[i] * log2(test_pro_ref[i]);
  }
  
  //std::cout<<"H_ref from b spline"<<H_ref_<<", H_ref from original definition "<<test_H_ref<<std::endl;
  
};

void EdgeSE3ProjectIntensityOnlyPoseNID::ClearPrevH(){
  pro_current_ = std::vector<float>(bin_num_,0.0);
  pro_ref_ = std::vector<float>(bin_num_,0.0);
  pro_joint_ = std::vector<std::vector<float> >(bin_num_);
  for(int i = 0; i < bin_num_; i++){
    pro_joint_[i] = std::vector<float>(bin_num_, 0.0);
  }
  H_joint_ = 0.0;
  H_current_ = 0.0;  
};

//pose update only for intensity

// bool EdgeSE3ProjectIntensityOnlyPose::read(std::istream& is){
//   return true;
// }

// bool EdgeSE3ProjectIntensityOnlyPose::write(std::ostream& os) const {
//   return os.good();
// }

// void EdgeSE3ProjectIntensityOnlyPose::linearizeOplus(){
//   //2d pixel position in l frame
//   float x = fx_ * u_l_ + cx_;
//   float y = fy_ * v_l_ + cy_;

//   //interpolation, DSO, getInterpolatedElement33() function, bilinear interpolation
//   int ix = (int)x;
//   int iy = (int)y;
//   float dx = x - ix;
//   float dy = y - iy;
//   float dxdy = dx*dy;

//   float gradient_x_11, gradient_x_10, gradient_x_01, gradient_x_00;
//   float gradient_y_11, gradient_y_10, gradient_y_01, gradient_y_00;
//   float gradient_x, gradient_y;

//   //calculate the four corner pixel's gradient and do bilinear inerpolation
//   if(ix >= 0 && ix+2 <= image_last_.cols - 1 && iy >= 0 && iy+2 <= image_last_.rows - 1){
//     gradient_x_11 = image_last_.ptr<float>(iy+1)[ix+2] - image_last_.ptr<float>(iy+1)[ix+1];
//     gradient_y_11 = image_last_.ptr<float>(iy+2)[ix+1] - image_last_.ptr<float>(iy+1)[ix+1];

//     gradient_x_10 = image_last_.ptr<float>(iy+1)[ix+1] - image_last_.ptr<float>(iy+1)[ix];
//     gradient_y_10 = image_last_.ptr<float>(iy+2)[ix] - image_last_.ptr<float>(iy+1)[ix];

//     gradient_x_01 = image_last_.ptr<float>(iy)[ix+2] - image_last_.ptr<float>(iy)[ix+1];
//     gradient_y_01 = image_last_.ptr<float>(iy+1)[ix+1] - image_last_.ptr<float>(iy)[ix+1];

//     gradient_x_00 = image_last_.ptr<float>(iy)[ix+1] - image_last_.ptr<float>(iy)[ix];
//     gradient_y_00 = image_last_.ptr<float>(iy+1)[ix] - image_last_.ptr<float>(iy)[ix];

//     gradient_x = dxdy * gradient_x_11 + (dy - dxdy) * gradient_x_10 + (dx - dxdy) * gradient_x_01 + (1 - dx - dy + dxdy) * gradient_x_00;
//     gradient_y = dxdy * gradient_y_11 + (dy - dxdy) * gradient_y_10 + (dx - dxdy) * gradient_y_01 + (1 - dx - dy + dxdy) * gradient_y_00;
//     //std::cout<<"gradient x and y "<<gradient_x<<", "<<gradient_y<<std::endl;
//   }
//   else{
//     //may need think more about this
//     //std::cout<<"opps..... out of image boundary"<<std::endl;
//     gradient_x = 0.0;
//     gradient_y = 0.0;
//     this->setLevel(1);
//   }

//   _jacobianOplusXi(0,0) = gradient_x * inv_depth_ * fx_;
//   _jacobianOplusXi(0,1) = gradient_y * inv_depth_ * fy_;
//   _jacobianOplusXi(0,2) = -inv_depth_ * (gradient_x * fx_ * u_l_ - gradient_y * fy_ * v_l_);
//   _jacobianOplusXi(0,3) = -gradient_x * fx_ * u_l_ * v_l_ - gradient_y * fy_ * (1.0 + v_l_*v_l_);
//   _jacobianOplusXi(0,4) = gradient_x * fx_ * (1.0 + u_l_ * u_l_) + gradient_y * fy_ * u_l_ * v_l_;
//   _jacobianOplusXi(0,5) = -gradient_x * fx_ * v_l_ + gradient_y * fy_ * u_l_;

//   //std::cout<<"the jacobian is \n"<<_jacobianOplusXi<<std::endl;
//   //_jacobianOplusXi = -_jacobianOplusXi;

//   if(isnan(_jacobianOplusXi(0,0)) || isnan(_jacobianOplusXi(0,1)) || isnan(_jacobianOplusXi(0,2)) || isnan(_jacobianOplusXi(0,3)) || isnan(_jacobianOplusXi(0,4)) || isnan(_jacobianOplusXi(0,5))){
//     this->setLevel(1);
//   }
// }


void EdgeSE3ProjectDirectOnlyPose::linearizeOplus(){
    if ( level() == 1 )
    {
        _jacobianOplusXi = Eigen::Matrix<double, 1, 6>::Zero();
        return;
    }
    VertexSE3Expmap* vtx = static_cast<VertexSE3Expmap*> ( _vertices[0] );
    Eigen::Vector3d xyz_trans = vtx->estimate().map ( x_world_ );   // q in book

    double x = xyz_trans[0];
    double y = xyz_trans[1];
    double invz = 1.0/xyz_trans[2];
    double invz_2 = invz*invz;

    float u = x*fx_*invz + cx_;
    float v = y*fy_*invz + cy_;

    // jacobian from se3 to u,v
    // NOTE that in g2o the Lie algebra is (\omega, \epsilon), where \omega is so(3) and \epsilon the translation
    Eigen::Matrix<double, 2, 6> jacobian_uv_ksai;

    jacobian_uv_ksai ( 0,0 ) = - x*y*invz_2 *fx_;
    jacobian_uv_ksai ( 0,1 ) = ( 1+ ( x*x*invz_2 ) ) *fx_;
    jacobian_uv_ksai ( 0,2 ) = - y*invz *fx_;
    jacobian_uv_ksai ( 0,3 ) = invz *fx_;
    jacobian_uv_ksai ( 0,4 ) = 0;
    jacobian_uv_ksai ( 0,5 ) = -x*invz_2 *fx_;

    jacobian_uv_ksai ( 1,0 ) = - ( 1+y*y*invz_2 ) *fy_;
    jacobian_uv_ksai ( 1,1 ) = x*y*invz_2 *fy_;
    jacobian_uv_ksai ( 1,2 ) = x*invz *fy_;
    jacobian_uv_ksai ( 1,3 ) = 0;
    jacobian_uv_ksai ( 1,4 ) = invz *fy_;
    jacobian_uv_ksai ( 1,5 ) = -y*invz_2 *fy_;

    Eigen::Matrix<double, 1, 2> jacobian_pixel_uv;

    jacobian_pixel_uv ( 0,0 ) = ( get_interpolated_pixel_value ( u+1,v )-get_interpolated_pixel_value ( u-1,v ) ) /2;
    jacobian_pixel_uv ( 0,1 ) = ( get_interpolated_pixel_value ( u,v+1 )-get_interpolated_pixel_value ( u,v-1 ) ) /2;

    _jacobianOplusXi = jacobian_pixel_uv*jacobian_uv_ksai;  
}

} // end namespace
