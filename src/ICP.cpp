#include <Eigen/Dense>
#include "slamBase.h"
typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Vertices;
Eigen::Matrix<double,4,4> Rt;
using namespace Eigen;
Vertices frame2vertices(FRAME f)
{
   Vertices vertices;
    for (int m = 0; m < f.depth.rows; m+=2)
        for (int n=0; n < f.depth.cols; n+=2)
        {
            // 获取深度图中(m,n)处的值
            ushort d = depth.ptr<ushort>(m)[n];
            // d 可能没有值，若如此，跳过此点
            if (d == 0)
                continue;
            // d 存在值，则向点云增加一个点
            PointT p;

            // 计算这个点的空间坐标
            p.z = double(d) / camera.scale;
            p.x = (n - camera.cx) * p.z / camera.fx;
            p.y = (m - camera.cy) * p.z / camera.fy;
            
           // 从rgb图像中获取它的颜色
            // rgb是三通道的BGR格式图，所以按下面的顺序获取颜色
/*          p.b = rgb.ptr<uchar>(m)[n*3];
            p.g = rgb.ptr<uchar>(m)[n*3+1];
            p.r = rgb.ptr<uchar>(m)[n*3+2];  
*/ 
	    //将三维坐标提到vertices
	    int curr_vertex = 0;
	    vertices(0,curr_vertex) = p.x;
	    vertices(1,curr_vertex) = p.y;
	    vertices(2,curr_vertex) = p.z;
	    curr_vertex++;
	}
	return vertices;

}
//homogeneous3d->2d
template<typename t,Index i,Index n>
Matrix<t,i-1,n> dehomogenize(Matrix<t,i,n> source)
{
  Matrix<t,2,n> output;
  for (Index j = 0;j< i-1;++j) output(j) = source(j)/source(i-1) ;
  return output;
}
//矩阵置为0或1
template<typename t,Index n,Index m>
Matrix<t,m,n> init0_1(Matrix<t,m,n> Mat,int it)
{
  for(Index j = 0; j < m; ++j)
     for(Index i = 0; i < n; ++i)
       Mat(j,i) = it;
     return Mat;
}
//get int image coord, then get its depth 
Vector3d find_assoc_3Dpoint(FRAME f,Vector2d id,CAMERA_INTRINSIC_PARAMETERS camera)
{
  Vector3d output;
  Vector2i iid;
  iid.x() = round(id.x());
  iid.y() = round(id.y());
  if(iid.x() >= f.depth.cols() || iid.y() >= f.depth.rows())
	return init0_1(output,0);
  double d = f.depth.Ptr<double>(iid.x())[iid.y()];//能否取到？iid超出边界如何处理？

  output(2) = d / camera.scale;
  output(0) = (f.depth.cols() - camera.cx) * output(2) / camera.fx;
  output(1) = (f.depth.rows() - camera.cy) * output(2) / camera.fy;
 // output(3) = 1;
  return output;
}
Eigen::Matrix<double,4,4> motionEstimate(FRAME frame0,FRAME frame1,Transform transform2model,CAMERA_INTRINSIC_PARAMETERS camera)
{
  //世界坐标
  Vertices Vset0 = frame2vertices(frame0);
  Vertices Vset1 = frame2vertices(frame1);
  //E = sum ((pi-qi)*ni+r*(pixni)+t*ni)^2;
  double Errors = 0;
  //投影法找匹配点
  //×
  for (int i = 0; i < Vset0.size;i++)
  {
    //计算由Vset0到1的transform
    Eigen::Matrix<double,2,1> model_idx= dehomogenize(transform2model(Vset0(i)));
    //model_p 来源于Vset1
    Vector3d model_p= find_assoc_3Dpoint(frame1, model_idx,camera);
    //Vector3d model_n = find_assoc_3dpoint(model_normals, model_idx);
   // model_n = model_n.normalized();
    Vector3d deviation = Vset0(i)-model_p;
    
  }
    
}
