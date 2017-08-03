
#include <Eigen/Dense>
#include <Eigen/Core>
#include <boost/graph/graph_concepts.hpp>
#include "slamBase.h"
#include "ICP.h"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
//#define DEBUG // 启用 IcpDEBUG 宏时将开启调试输出
//#define DEBUG_ITER
#ifdef DEBUG_ITER
#include <stdio.h>
template <typename S,int R>
void IcpDEBUG_iter(Eigen::Matrix<S,R,1> V)
{
  std::cout<<V<<'\n';
}
template <typename S>
void IcpDEBUG_iter(S a)
{
  std::cout<<a<<'\n';
}
#else
#define IcpDEBUG_iter(b)
#endif


#ifdef  DEBUG
#include <stdio.h>
template <typename S>
void IcpDEBUG(S a)
{
  std::cout<<a<<'\n';
}
void IcpDEBUG()
{
  printf("调试输出结束位..\n");
}
#else
#define  IcpDEBUG(a)
#endif

typedef Eigen::Matrix<double,3,1> Vector3d;
typedef Eigen::Matrix<Vertice, Eigen::Dynamic,1> Vertices;//还需测试
using namespace Eigen;

void BilateralFilter_depth(cv::Mat depth,cv::Mat output)
{
  const int N = 4;
  const double theta_r_1 = 0.0285714;//34
  const double theta_s_1 = 0.2;//5
  for(int m = 4;m < depth.rows; m++)
    for(int n = 4;n < depth.cols; n++)
    {
      double w_s[m+N][n+N],w_r[m+N][n+N],w[m+N][n+N];
      double w_total = 0;
      double wxg_total = 0;
      for(int i = m-N; i < m+N;i++)
	for (int j = n-N; j < n+N;j++)
	{
	  if(i < 0||j < 0||i >= depth.rows || j >= depth.cols)
	    continue;
	  w_s[i][j] = exp(((i-m)*(i-m)+(j-n)*(j-n))*0.5*theta_s_1*theta_s_1);
	  double diff = double(depth.ptr<ushort>(i)[j] - depth.ptr<ushort>(m)[n])*0.05;//先缩小20倍，否则w_r可能会超界，e3000*3000
	  w_r[i][j] = exp(diff*diff*0.5*theta_r_1*theta_r_1);
	  w[i][j] = w_s[i][j]*w_r[i][j];
	  w_total += w[i][j];
	  wxg_total += w[i][j]*double(depth.ptr<ushort>(i)[j])*0.05;
	}
	double dp = wxg_total/w_total;
      output.ptr<double>(m)[n] = dp*20;
    }
}

//get int image coord, then get its depth 
//(u,v)对应(n,m)!!!!!!
Vector3d find_assoc_3Dpoint(int n,int m,cv::Mat depth,CAMERA_INS camera)//越界返回的是0
{
    if(m >= depth.rows || n >= depth.cols || m < 1 || n < 1)
       return Vector3d::Zero();
    auto d = depth.ptr<ushort>(m)[n];//注意类型，
    if(d == 0)
       return Vector3d::Zero();
    // 计算这个点的空间坐标,1/fx = 0.0019305,1/fy = 0.00192678,(u,v),u=n,v=m
    double z = double(d) *0.001;
    double x = (double(n) - camera.cx) * z * 0.0019305;
    double y = (double(m) - camera.cy) * z * 0.00192678;  
    return Vector3d(x,y,z);
}

Vector3d find_assoc_3Dpoint(Vector2i iid,cv::Mat depth,CAMERA_INS camera)//越界返回的是0
{
 return find_assoc_3Dpoint(iid(0),iid(1), depth, camera);
}
//在映射点区域搜索最邻近点
Vector3d find_closest_assoc_3Dpoint(Vector3d target,Vector2i imgid,cv::Mat target_depth,CAMERA_INS camera)
{
  const int ksize = 1;//在5x5的kernel里找
  int r = imgid(0);
  int c = imgid(1);
  int sized = (2*ksize+1)*(2*ksize+1);
  double dist[sized];
  int num = 0;
  //Vector3d certer = find_assoc_3Dpoint(imgid,target_depth,camera);
  for (int dr = r-ksize; dr <= r+ksize; dr++)
    for(int dc = c-ksize; dc <= c+ksize; dc++)
    {
      Vector3d dpoint = find_assoc_3Dpoint(dr,dc,target_depth,camera);
      dist[num] = (target-dpoint).norm();
      num++;
    }
  int position = min_element(dist,dist+sized) - dist;//位置差
  int m = r-ksize + position/(2*ksize+1);
  int n = c-ksize + position - (m-r+ksize)*(2*ksize+1);
  return find_assoc_3Dpoint(m,n,target_depth,camera);
}
//深度差分得normal
Vector3d compute_m_normal(Vector2i iid,cv::Mat depth,CAMERA_INS camera)
{
    const int kHalfRegion = 2;
    Vector2i tmp1,tmp2;
    Vector3d dx = Vector3d::Zero(), dy = Vector3d::Zero();
    for (int y = -kHalfRegion; y <= kHalfRegion; ++y) 
    {
      tmp1 = iid + Vector2i(-kHalfRegion, y);
      tmp2 = iid + Vector2i(kHalfRegion, y);
      Vector3d p1 = find_assoc_3Dpoint(tmp1,depth,camera);
      Vector3d p2 = find_assoc_3Dpoint(tmp2,depth,camera);
      //有一个出界的话，normal就不准，都不要
      if(p1 == Vector3d::Zero() || p2 == Vector3d::Zero())
      {
	 return Vector3d::Zero();//没有normal情况返回0，上层处理
	 break;
      }
      dx += Vector3d(p2(0) - p1(0), p2(1) - p1(1), p2(2) - p1(2));//
      }
    for (int x = -kHalfRegion; x <= kHalfRegion; ++x) 
    {
      	  tmp1 = iid + Vector2i(x, -kHalfRegion);
	  tmp2 = iid + Vector2i(x, kHalfRegion);
	  Vector3d p1 = find_assoc_3Dpoint(tmp1,depth,camera);
	  Vector3d p2 = find_assoc_3Dpoint(tmp2,depth,camera);
	  //有一个出界的话，normal就不准，都不要
	  if(p1 == Vector3d::Zero() || p2 == Vector3d::Zero())
	  {
		return Vector3d::Zero();//没有normal情况返回0，上层处理
		break;
	  }
	  dy += Vector3d(p2(0) - p1(0), p2(1) - p1(1), p2(2) - p1(2));
       }  
    Vector3d output = (dx.cross(dy)).normalized();
    if (output(2) > 0.f) output = -output;
    return output;
   // normal[idx] = Vector4f(n.x(), n.y(), n.z(), 1);
}
//3d point project->2d image coord
Vector2d project2(Matrix<double,3,1> point3,CAMERA_INS camera)
{
  //x = fx*X + Cx*Z ; y = fy*Y + Cy*Z ; z = Z; 
  Vector2d output; //(u,v)
  output(0) = camera.fx*point3(0)/point3(2) + camera.cx;
  output(1) = camera.fy*point3(1)/point3(2) + camera.cy;
//  IcpDEBUG_iter(output);
  return output;
}
bool useless(Vector3d p)
{
  //差不多以图像长宽边界，然后太近太远的点也不要
  if(fabs(p(0) > 3) || fabs(p(0)) < 0.01|| fabs(p(1) > 2) ||fabs(p(1)) < 0.01|| p(2) > 10 ||p(2) < 0.2) 
    return true;
  else
    return false;
}

RigidTransf<double> compute_rigidTransform (RigidTransf<double> transform2model,double & E_linear,FRAME f0,FRAME f1,CAMERA_INS camera,double thresh_dis,double thresh_angle)
{
    RigidTransf<double>  outZero;
    Eigen::Matrix<double,6,6> ATA = Eigen::Matrix<double,6,6>::Zero();
    Eigen::Matrix<double,6,1> ATb = Eigen::Matrix<double,6,1>::Zero();
    Eigen::Matrix<double,6,1> cof_x;
    int pair_index = 0;
    for (int i = 0; i < f0.depth.rows; i++)//u
    for (int j = 0; j < f0.depth.cols; j++)//v
      {
	//计算由Vset0到1的transform
	Vector3d local_p0 = find_assoc_3Dpoint(j,i,f0.depth,camera);//(u,v)
	Vector3d local_p = transform2model.Transf(local_p0);
	if(useless(local_p))
	   continue;
	//这里近似的把点在两个图像中的坐标看做一样的，model_idx即是frame0,也是frame1中，而model_p来源于Vset1
	Vector2d model_ixd= project2(local_p,camera);//返回的可能是一个越界的图像坐标，下面得到的p,n就会都是0了。
	Vector2i model_ix = Vector2i(round(model_ixd(0)),round(model_ixd(1)));
	//取整后，去映射点小区域搜索最邻近点
	Vector3d model_p= find_closest_assoc_3Dpoint(local_p,model_ix,f1.depth,camera);
	Vector3d model_n = compute_m_normal(model_ix,f1.depth,camera);
	if(useless(model_p)|| model_n == Vector3d::Zero())  //直接跳过，好像会引导优化方向都跳过，会不收敛
	    continue;
	//IcpDEBUG_iter(model_n);
	//距离上限约束和方向一致性约束
	double distances = (local_p - model_p).norm();
	double angle_cos = (local_p.dot(model_n))/local_p.norm();//angle 控制在30度之内，angle_cos > cos 30
	//IcpDEBUG_iter(distances);
	IcpDEBUG_iter(angle_cos);
	if(distances > thresh_dis ||fabs(angle_cos) < thresh_angle)
	   continue;
    /* 
     * 线性化，首先pi,qi坐标齐次化,pi==Vset0i, qi==model_p,
     * 
     * 
     E = sum [(M*pi -qi)*ni]^2
     
	      [1  -y  B  tx]
     M = T*R =|y   1 -a  ty|
	      |-B  a  1  tz|
	      [0   0  0   1]
    令x = （a, B, y, tx, ty, tz）T
    E = min |Ax -b|^2
    A = n'*G , G = [Px | I3 ],Px 为pi的坐标反对称矩阵，I3为3x3单位矩阵，n 为法向量,A最后为1x6矩阵,M个匹配点对就是Mx6
    =>A = [nz*Py-ny*Pz, nx*Pz-nz*Px, ny*Px-nx*Py, nx, ny, nz]
    b = n' *(P-Q)
    =>b = [nx*Qx+ny*Qy+nz*Qz-nx*Px-ny*Py-nz*Pz] 

	Eigen::Matrix<double,1,6> A_i = Eigen::Matrix<double,1,6>::Identity();
	A_i(0) = model_n(2)*local_p(1) - model_n(1)*local_p(2); //0=x;1=y;2=z;
	A_i(1) = model_n(0)*local_p(2) - model_n(2)*local_p(0);
	A_i(2) = model_n(1)*local_p(0) - model_n(0)*local_p(1);
	A_i(3) = model_n(0);
	A_i(4) = model_n(1);
	A_i(5) = model_n(2);
	double b_i = model_n(0)*( model_p(0)-local_p(0) ) + model_n(1)*( model_p(1)-local_p(1) ) + model_n(2)*( model_p(2)-local_p(2) );
// 	if(fabs(b_i) < 0.0001 || b_i > 1e+100 ||std::isnan(b_i))
// 	   continue;
	cof_A.conservativeResize(pair_index + 1, 6);//riSize出了问题，它把之前的值都改变了
	cof_b.conservativeResize(pair_index + 1);
	for(int i = 0; i < 6;i++)
	    cof_A(pair_index,i) = A_i(i);
	cof_b(pair_index) = b_i;
	++pair_index;
  以上是之前采用eigen ls solver 的方法。
 * 
 */
	const float  sx = local_p(0);
	const float  sy = local_p(1);
	const float  sz = local_p(2);
	const float  dx = model_p(0);
	const float  dy = model_p(1);
	const float  dz = model_p(2);
	const float  nx = model_n(0);
	const float  ny = model_n(1);
	const float  nz = model_n(2);
	double a = nz*sy - ny*sz;
	double b = nx*sz - nz*sx; 
	double c = ny*sx - nx*sy;
    //    0  1  2  3  4  5
    //    6  7  8  9 10 11
    //   12 13 14 15 16 17
    //   18 19 20 21 22 23
    //   24 25 26 27 28 29
    //   30 31 32 33 34 35
	ATA  (0,0) += a * a;
	ATA  (0,1) += a * b;
	ATA  (0,2) += a * c;
	ATA  (0,3) += a * nx;
	ATA  (0,4) += a * ny;
	ATA  (0,5) += a * nz;
	ATA  (1,1) += b * b;
	ATA  (1,2) += b * c;
	ATA  (1,3) += b * nx;
	ATA  (1,4) += b * ny;
	ATA  (1,5) += b * nz;
	ATA  (2,2) += c * c;
	ATA  (2,3) += c * nx;
	ATA  (2,4) += c * ny;
	ATA  (2,5) += c * nz;
	ATA  (3,3) += nx * nx;
	ATA  (3,4) += nx * ny;
	ATA  (3,5) += nx * nz;
	ATA  (4,4) += ny * ny;
	ATA  (4,5) += ny * nz;
	ATA  (5,5) += nz * nz;
        double d = nx*dx + ny*dy + nz*dz - nx*sx - ny*sy - nz*sz;
	ATb(0) += a * d;
	ATb(1) += b * d;
	ATb(2) += c * d;
	ATb(3) += nx * d;
	ATb(4) += ny * d;
	ATb(5) += nz * d;
	pair_index++;
	Eigen::Matrix<double,3,2> lm;
	lm.col(0) = local_p;
	lm.col(1) = model_p;
 	IcpDEBUG_iter(lm);
//为了方便调试，简化	
//  	if(pair_index > 1000)
//  	    break;
      }
      IcpDEBUG();
      IcpDEBUG(pair_index);
//    线性最小二乘问题，normal equations 最快，QR次之，SVD最慢。
//    问题a，B，y都是欧拉角，需要加优化的约束吗？
//    Eigen::Matrix<double,6,1> cof_x = (cof_A.transpose()*cof_A).ldlt().solve(At_b);
//    IcpDEBUG(cof_A);
//    IcpDEBUG();
      if(pair_index < 20)//点太少 可能无解
	return outZero;
     cof_x = (ATA.inverse() * ATb);
     Matrix3d R = eulerAngle2R(cof_x(0),cof_x(1),cof_x(2));
     Vector3d T = Vector3d(cof_x(3),cof_x(4),cof_x(5));
     //误差
     E_linear = (ATA*cof_x - ATb).norm()/ATb.norm();
     IcpDEBUG(E_linear);    
     IcpDEBUG(Vector3d(cof_x(0),cof_x(1),cof_x(2)));
     IcpDEBUG(T);
     RigidTransf<double>  Trans2Model(R,T);
    // std::cout<<"中间变换矩阵transform"<<Trans2Model<<std::endl;
     return Trans2Model;
  
}
Isometry3d motionEstimate(FRAME frame0, FRAME frame1,CAMERA_INS camera, int icp_iter,double thresh_dis,double thresh_angle)
{
  //世界坐标,需不需要双边滤波？暂时不
 // Vertices Vset0 = frame2vertices(frame0,camera);
 // Vertices Vset1 = frame2vertices(frame1,camera);
//双边滤波
//   cv::Mat depthtmp0(frame0.depth.size(),6);
//   BilateralFilter_depth(frame0.depth,depthtmp0);
//   frame0.depth = depthtmp0;
//   cv::Mat depthtmp1(frame1.depth.size(),6);
//   BilateralFilter_depth(frame1.depth,depthtmp1);
//   frame1.depth = depthtmp1;
  int m = frame0.depth.rows;
  int n = frame0.depth.cols;
  if(m*n < 100)//点少的图像不能做估计,边界上的点也不要
       return Isometry3d::Identity();
  //原始的为非线性 E = sum[(Rpi+t-qi)*ni]^2
  //线性化： E = sum [(pi-qi)*ni+r*(pixni)+t*ni]^2;
  //double E_nonlinear = 0;
  Matrix<double,20,1> E_list;
  RigidTransf<double> trans2model;
  //ICP 迭代次数
  for(int iter = 0; iter < icp_iter; iter++)
  {
    double E_linear = 0;
    RigidTransf<double> trans_incre = compute_rigidTransform(trans2model,E_linear,frame0,frame1,camera,thresh_dis,thresh_angle);
    trans2model.R = trans2model.R * trans_incre.R;
    trans2model.t = trans2model.t + trans_incre.t;
    E_list(iter) = E_linear;
  }
  Isometry3d _transform = Isometry3d::Identity();
  _transform.rotate(trans2model.R);
  _transform.translation() = trans2model.t;
   return _transform;
}
