
#include <Eigen/Dense>
#include <Eigen/Core>
#include "slamBase.h"
#include "ICP.h"

#ifndef __TESTDBG_H__
#define __TESTDBG_H__

#define DBG // 启用 DBG 宏时将开启调试输出

#ifdef  DBG
#include <stdio.h>
template <typename S>
void testDBG(S a)
{
  std::cout<<a<<'\n';
}
template <typename S,int R>
void testDBG(Eigen::Matrix<S,R,1> V)
{
  std::cout<<V<<'\n';
}

void testDBG()
{
  printf("调试输出结束位..\n");
}
#else
#define  testDBG(a)
#endif
#endif

typedef Eigen::Matrix<double,3,1> Vector3d;
typedef Eigen::Matrix<Vertice, Eigen::Dynamic,1> Vertices;//还需测试
using namespace Eigen;

double thresh_dis = 2.7;
Vertices frame2vertices(FRAME f,CAMERA_INS camera)//
{
   Vertices Vsets;
   int num = 0;
    for (int m = 0; m < f.depth.rows; m+=2)
    {
      for (int n = 0; n < f.depth.cols; n+=2)
        {
            // 获取深度图中(m,n)处的值,m为行，n为列
            ushort d = f.depth.ptr<ushort>(m)[n];
            // d 可能没有值，若如此，跳过此点
            if (d == 0)
                continue;
            // d 存在值，则向点云增加一个点
            Vector3d p;
            // 计算这个点的空间坐标
            p(2) = double(d) / camera.scale;
            p(0) = (n - camera.cx) * p(2) / camera.fx;
            p(1) = (m - camera.cy) * p(2) / camera.fy;  
	   // testDBG(p);
           // 从rgb图像中获取它的颜色
            // rgb是三通道的BGR格式图，所以按下面的顺序获取颜色
/*          p.b = rgb.ptr<uchar>(m)[n*3];
            p.g = rgb.ptr<uchar>(m)[n*3+1];
            p.r = rgb.ptr<uchar>(m)[n*3+2];  
*/ 
	    //将三维坐标output提到vertices。对应的图像坐标放到最后两位。
	    Vsets.resize(num+1);
	    Vsets(num).worldCoord = p;
	    Vsets(num).imgCoord = Vector2i(m,n);
	    ++num;
	}
    }
    return Vsets;

}

//get int image coord, then get its depth 
Vector3d find_assoc_3Dpoint(int m,int n,Vertices Vset,FRAME f,CAMERA_INS camera)//越界返回的是0
{
    Vector3d p;
    if(m >= f.depth.rows || n >= f.depth.cols)
       return Vector3d::Zero();
    ushort d = f.depth.ptr<ushort>(m)[n];
    if(d == 0)
       return Vector3d::Zero();
    // 计算这个点的空间坐标
    p(2) = double(d) / camera.scale;
    p(0) = (double(n) - camera.cx) * p(2) / camera.fx;
    p(1) = (double(m) - camera.cy) * p(2) / camera.fy;  
    return p;
}
Vector3d find_assoc_3Dpoint(Vector2i iid,Vertices Vset,FRAME f,CAMERA_INS camera)//越界返回的是0
{
  Vector3d out = find_assoc_3Dpoint(iid(0),iid(1),Vset, f, camera);
  return out;
}
//深度差分得normal
Vector3d compute_m_normal(Vector2i iid,Vertices Vset,FRAME f,CAMERA_INS camera)
{
    const int kHalfRegion = 2;
    Vector2i tmp1,tmp2;
    Vector3d dx = Vector3d::Zero(), dy = Vector3d::Zero();
    for (int y = -kHalfRegion; y <= kHalfRegion; ++y) 
    {
      tmp1 = iid + Vector2i(-kHalfRegion, y);
      tmp2 = iid + Vector2i(kHalfRegion, y);
      Vector3d p1 = find_assoc_3Dpoint(tmp1,Vset,f,camera);
      Vector3d p2 = find_assoc_3Dpoint(tmp2,Vset,f,camera);
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
	  Vector3d p1 = find_assoc_3Dpoint(tmp1,Vset,f,camera);
	  Vector3d p2 = find_assoc_3Dpoint(tmp2,Vset,f,camera);
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
Matrix<double,2,1> project2(Matrix<double,3,1> point3,CAMERA_INS camera)
{
  //x = fx*X + Cx*Z ; y = fy*Y + Cy*Z ; z = Z;
  Matrix<double,2,1> output;
  output(0) = camera.fx*point3(0)/point3(2) + camera.cx;
  output(1) = camera.fy*point3(1)/point3(2) + camera.cy;
  return output;
}
bool useless(Vector3d p)
{
  if(fabs(p(0) > 3) || fabs(p(0)) < 0.01|| fabs(p(1) > 1) ||fabs(p(1)) < 0.01|| p(2) > 3 ||p(2) < 0.2) //太近太远的点都不要
    return true;
  else
    return false;
}
RigidTransf<double> compute_rigidTransform (Vertices Vset0, Vertices Vset1,RigidTransf<double> transform2model,double & E_linear,FRAME f1,CAMERA_INS camera)
{
    Matrix<PointPair,Eigen::Dynamic,1> pairs;//储存匹配点信息
    Matrix<double, Eigen::Dynamic,6> cof_A;//需要自己回收吗？
    Matrix<double, Eigen::Dynamic,1> cof_b;
    int pair_index = 0;
    for (int i = 0; i < Vset0.size();i++)
      {
	//计算由Vset0到1的transform
	Vector3d local_p = Vset0(i).worldCoord;
	if(useless(local_p))
	   continue;
	//testDBG(Vset0(i).worldCoord);
	Vector2d model_ixd= project2(transform2model.Transf(local_p),camera);//返回的可能是一个越界的图像坐标，下面得到的p,n就会都是0了。
	Vector2i model_ix = Vector2i(round(model_ixd(0)),round(model_ixd(1)));
	//这里近似的把点在两个图像中的坐标看做一样的，model_idx即是frame0,也是frame1中，而model_p来源于Vset1
	Vector3d model_p= find_assoc_3Dpoint(model_ix,Vset1,f1,camera);
	Vector3d model_n = compute_m_normal(model_ix,Vset1,f1,camera);
	if(useless(model_p) || model_n == Vector3d::Zero())  //怎么处理？？？,直接跳过，好像会引导优化方向都跳过，会不收敛
	    continue;
	//testDBG(model_n);
	//距离上限约束和方向一致性约束
	double distances = (local_p - model_p).norm();
	double angle_cos = (local_p.dot(model_n))/local_p.norm();//angle 控制在60度之内，angle_cos > cos 60
	//testDBG(distances);
	//testDBG(angle_cos);
	if(distances > thresh_dis ||fabs(angle_cos) < 0.5)
	   continue;
//   E_nonlinear +=  Dot(deviation,model_n); 
    /* 线性化，首先pi,qi坐标齐次化,pi==Vset0i, qi==model_p,
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
    */
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
	cof_A.conservativeResize(pair_index + 1, 6);//出了问题，好像他把之前的值都改变了
	cof_b.conservativeResize(pair_index + 1);
	for(int i = 0; i < 6;i++)
	    cof_A(pair_index,i) = A_i(i);
	cof_b(pair_index) = b_i;
	PointPair pa(local_p,model_p,model_ix,model_n);//依次是local，model，idx，model_n
	pairs.resize(pair_index+1,1);
	pairs(pair_index) = pa;
	++pair_index;
	
// 	testDBG(local_p);
// 	testDBG(model_p);
// 	testDBG(model_n);
//为了方便调试，简化	
//  	if(pair_index > 1000)
//  	    break;
      }
      std::cout<<"匹配点对数:"<<pair_index<<std::endl;
// 线性最小二乘问题，normal equations 最快，QR次之，SVD最慢。
// 问题a，B，y都是欧拉角，需要加优化的约束吗？
//      Eigen::Matrix<double,6,1> cof_x = (cof_A.transpose()*cof_A).ldlt().solve(At_b);
 //    testDBG(cof_A);
 //    testDBG();
     Eigen::Matrix<double,6,1> cof_x = cof_A.colPivHouseholderQr().solve(cof_b);
     Matrix3d R = eulerAngle2R(cof_x(0),cof_x(1),cof_x(2));
     Vector3d T = Vector3d(cof_x(3),cof_x(4),cof_x(5));
//      if(std::isnan(R(0,0)) || std::isnan(T(0)))
//        return RigidTransf<double> zero;
     testDBG(Vector3d(cof_x(0),cof_x(1),cof_x(2)));
     testDBG(T);
    // testDBG(cof_b);
     RigidTransf<double>  Trans2Model(R,T);
     //误差
     E_linear = (cof_A*cof_x - cof_b).norm();
     std::cout<<"中间误差："<<E_linear<<'\n';
    // std::cout<<"中间变换矩阵transform"<<Trans2Model<<std::endl;
     return Trans2Model;
  
}
Isometry3d motionEstimate(FRAME frame0, FRAME frame1,CAMERA_INS camera)
{
  //世界坐标,需不需要双边滤波？暂时不
  Vertices Vset0 = frame2vertices(frame0,camera);
  Vertices Vset1 = frame2vertices(frame1,camera);
  if(Vset0.size() < 20 || Vset1.size() < 20)//点少的图像不能做估计,边界上的点也不要
       return Isometry3d::Identity();
  //原始的为非线性 E = sum[(Rpi+t-qi)*ni]^2
  //线性化： E = sum [(pi-qi)*ni+r*(pixni)+t*ni]^2;
  //double E_nonlinear = 0;
  Matrix<double,10,1> E_list;
  RigidTransf<double> trans2model;
  //ICP 迭代次数
  for(int iter = 0; iter < 10; iter++)
  {
    double E_linear = 0;
    RigidTransf<double> trans_incre = compute_rigidTransform(Vset0,Vset1,trans2model,E_linear,frame1,camera);
    trans2model.R = trans2model.R * trans_incre.R;
    //testDBG(trans_incre.t);
    trans2model.t = trans2model.t + trans_incre.t;
    E_list(iter) = E_linear;
  }
  Isometry3d _transform = Isometry3d::Identity();
  _transform.rotate(trans2model.R);
  _transform.translation() = trans2model.t;
   return _transform;
}
