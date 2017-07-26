#include <Eigen/Dense>
#include "slamBase.h"
#include "ICP.h"
typedef Eigen::Matrix<double,3,1> Vector3d;
typedef Eigen::Matrix<Vector3d, 2, Eigen::Dynamic> Vertices;//
//Eigen::Matrix<double,4,4> Rt;
using namespace Eigen;
Vertices frame2vertices(FRAME f)
{
   Vertices vsets;
    for (int m = 0; m < f.depth.rows; m+=1)
        for (int n=0; n < f.depth.cols; n+=1)
        {
            // 获取深度图中(m,n)处的值,m为行，n为列
            ushort d = depth.ptr<ushort>(m)[n];
            // d 可能没有值，若如此，跳过此点
            if (d == 0)
                continue;
            // d 存在值，则向点云增加一个点
            Vector3d p;

            // 计算这个点的空间坐标
            p(2) = double(d) / camera.scale;
            p(0) = (n - camera.cx) * p.z / camera.fx;
            p(1) = (m - camera.cy) * p.z / camera.fy;
            
           // 从rgb图像中获取它的颜色
            // rgb是三通道的BGR格式图，所以按下面的顺序获取颜色
/*          p.b = rgb.ptr<uchar>(m)[n*3];
            p.g = rgb.ptr<uchar>(m)[n*3+1];
            p.r = rgb.ptr<uchar>(m)[n*3+2];  
*/ 
	    //将三维坐标提到vertices。对应的图像坐标放到最后两位。
	    vsets(m,n) = p;	    
	}
	return vsets;

}

//get int image coord, then get its depth 
Vector3d find_assoc_3Dpoint(FRAME f,Vector2d id,CAMERA_INTRINSIC_PARAMETERS camera)
{
  Vector3d output;
  Vector2i iid(round(idx(0)),round(idx(1)));
  if(iid(0) >= f.depth.cols() || iid(1) >= f.depth.rows())
	return init_Mat(output,0);
  double d = f.depth.Ptr<double>(iid(0))[iid(1)];//能否取到？iid超出边界如何处理？
  output(2) = d / camera.scale;
  output(0) = (f.depth.cols() - camera.cx) * output(2) / camera.fx;
  output(1) = (f.depth.rows() - camera.cy) * output(2) / camera.fy;
 // output(3) = 1;
  return output;
}
//深度差分
Vector3d compute_m_normal(Vector2d idx,Vertices Vset,FRAME f)
{
    Vector3d output;
    const int kHalfRegion = 2;
    Vector2i iid(round(idx(0)),round(idx(1)));
    if(iid(0) >= f.depth.cols() || iid(1) >= f.depth.rows())
	return init_Mat(output,0);
    Vector3d dx = Vector3d::Zero(), dy = Vector3d::Zero();
    bool find_dx = false, find_dy = false;
    for (int y = -kHalfRegion; y <= kHalfRegion; ++y) 
    {
      if(iid(0)+kHalfRegion < m && iid(1)+y < n && iid(0)-kHalfRegion >= 0 && iid(1)+y >= 0) //越界则去掉,好像去掉的比他的多？
      {
        Vector3d p1 = Vset(iid + Vector2i(-kHalfRegion, y));
        Vector3d p2 = Vset(iid + Vector2i(kHalfRegion, y));
        //if (p1.w() == 0.f || p2.w() == 0.f) continue; 
        dx += Vector3d(p2(0) - p1(0), p2(1) - p1(1), p2(2) - p1(2));//不就等于【2,0.0】？
        find_dx = true;
      }
    }
    for (int x = -kHalfRegion; x <= kHalfRegion; ++x) 
    {
       if(iid(0)+x < m && iid(1)+kHalfRegion < n && iid(0)+x >= 0 && iid(1)-kHalfRegion >= 0)
       {
	  Vector3d p1 = Vset(iid + Vector2i(x, -kHalfRegion));
	  Vector3d p2 = Vset(iid + Vector2i(x, kHalfRegion));
	// if (p1.w() == 0.f || p2.w() == 0.f) continue;
	  dy += Vector3d(p2(0) - p1(0), p2(1) - p1(1), p2(2) - p1(2));
	  find_dy = true;
       }
    }

    if (!find_dx || !find_dy) return;//return ??? 没有normal情况应该跳过，具体处理？？

    Vector3d output = Cross(dx,dy).normalized();
    if (output(2) > 0.f) output = -output;
    return output;
   // normal[idx] = Vector4f(n.x(), n.y(), n.z(), 1);
}
Eigen::Matrix<double,4,4> motionEstimate(FRAME frame0,FRAME frame1,RigidTransf transform2model,CAMERA_INTRINSIC_PARAMETERS camera)
{
  //世界坐标,需不需要双边滤波？
  Vertices Vset0 = frame2vertices(frame0);
  Vertices Vset1 = frame2vertices(frame1);
  Eigen::Matrix<float, 6, 6> hessian;
  Eigen::Matrix<float, 6, 1> nabla;
  //原始的为非线性 E = sum[(Rpi+t-qi)*ni]^2
  //线性化： E = sum [(pi-qi)*ni+r*(pixni)+t*ni]^2;
  double E_nonlinear = 0;
  double E_linear = 0;
  for (int i = 0; i < Vset0.size();i++)
  {
    //计算由Vset0到1的transform
    Eigen::Matrix<double,2,1> model_idx= dehomogenize(transform2model.Transf(Vset0(i)));
    //这里近似的把点在两个图像中的坐标看做一样的，model_idx即是frame0,也是frame1中，而model_p来源于Vset1
    Vector3d model_p= find_assoc_3Dpoint(frame1, model_idx,camera);
    Vector3d model_n = compute_m_normal(model_idx,Vset1,frame1);
    Vector3d deviation =transform2model.Transf(Vset0(i)) - model_p;
    E_nonlinear +=  Dot(deviation,model_n); 
    
  //  double cos_angle = fabs(Dot(Vset0(i), model_n));
    
    
  }
    
}
