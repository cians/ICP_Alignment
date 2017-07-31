# pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
using namespace Eigen;
typedef int Index;
struct Vertice //点的世界坐标和图像坐标
{
  Vector3d worldCoord;
  Vector2i imgCoord;
};
//eulerAngle2RotationMatrix
Matrix3d eulerAngle2R(double alpha, double beta, double gamma) //alpha ->x; beta ->y; gamma ->z
{
  Matrix3d R;
  R(0,0) = 1;      R(0,1) = alpha*beta - gamma;  R(0,2) = alpha*gamma +beta;
  R(1,0) = gamma;  R(1,1) = alpha*beta*gamma+1;  R(1,2) = beta*gamma - alpha;
  R(2,0) = -beta;  R(2,1) = alpha;               R(2,2) =  1;
  return R;
}
//pairs
class PointPair
{
  public://默认初始化为0
    PointPair(Vector3d lo,Vector3d mo,Vector2i id,Vector3d mo_n):local(lo),model(mo),idx(id),modle_normal(mo_n){}
    PointPair(){}
    Vector3d local;
    Vector3d model;
    Vector2i idx;
    Vector3d modle_normal;
    
};
//刚体变换R,t
template <typename T>
class RigidTransf
{
  public:
    RigidTransf():R(Matrix<T,3,3>::Identity()),t(Matrix<T,3,1>::Zero()){}//R初始化为单位矩阵，t为0
    RigidTransf(Matrix<T,3,3> MR, Matrix<T,3,1> MT):R(MR),t(MT){}
    Matrix<T,3,1> Transf(Matrix<T,3,1> X) //将R,t变换应用到输入X上
    {
       return this->R*X + this->t;
    }
    Matrix<T,3,3> R;
    Matrix<T,3,1> t;
};
//反对称矩阵
template<typename S>
Matrix<S, 3, 3> skew_symmetry(Matrix<S,3,1> in)
{
  Matrix<S, 3, 3> out;
  out.row(0) = Vector3d(0, -in(2), in(1));
  out.row(1) = Vector3d(in(2), 0, -in(0));
  out.row(2) = Vector3d(-in(1), in(0), 0);
  return out;
}

//矩阵全置为it
template<typename t,Index n,Index m>
Matrix<t,m,n> init_Mat(Matrix<t,m,n> Mat,int it)
{
  for(Index j = 0; j < m; ++j)
     for(Index i = 0; i < n; ++i)
       Mat(j,i) = it;
     return Mat;
}

//矩阵相应位置相乘
template <typename S, int R>
S Dot( Matrix<S, R,1> &lhs, Matrix<S, R,1> &rhs)
{
    S result = 0;
    for (Index i = 0; i < R; ++i) result += lhs[i] * rhs[i];
    return result;
}