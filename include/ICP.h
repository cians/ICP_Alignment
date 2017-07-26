# pragma once
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
using namespace Eigen;
typedef int Index;
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

//homogeneous3d->2d
template<typename t,Index i,Index n>
Matrix<t,i-1,n> dehomogenize(Matrix<t,i,n> source)
{
  Matrix<t,2,n> output;
  for (Index j = 0;j< i-1;++j) output(j) = source(j)/source(i-1) ;
  return output;
}
//homogeneous
template<typename S,Index R>
Matrix<S, R + 1, 1> homogeneous() const
    {
        Matrix<S, R + 1, 1> h;
        for (Index r = 0; r < R; r++) h(r) = (*this)(r);
        h(R) = 1;
        return h;
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
//Cross
template <typename S>
_CPU_AND_GPU_ Vector<S, 3> Cross(const Vector<S, 3> &lhs,
                                 const Vector<S, 3> &rhs)
{
    Vector<S, 3> result;
    result.x() = lhs.y() * rhs.z() - lhs.z() * rhs.y();
    result.y() = -lhs.x() * rhs.z() + lhs.z() * rhs.x();
    result.z() = lhs.x() * rhs.y() - lhs.y() * rhs.x();
    return result;
}
//矩阵相应位置相乘
template <typename S, int R>
S Dot(const Vector<S, R> &lhs, const Vector<S, R> &rhs)
{
    S result = 0;
    for (Index i = 0; i < R; ++i) result += lhs[i] * rhs[i];
    return result;
}