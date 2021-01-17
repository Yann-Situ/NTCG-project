/* laplacian.cpp
 * Author : Yann-Situ GAZULL
 * Description :
 */
#include "laplacian.h"

#include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/eigs.h>
#include <iostream>

Lap::Lap(Eigen::MatrixXd _V, Eigen::MatrixXi _F)
: V(_V), F(_F)
{
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
    igl::cotmatrix(V,F,L);
}
Lap::~Lap(){}

Eigen::VectorXd Lap::curvature()
{
    Eigen::VectorXd K;
    // Compute integral of Gaussian curvature
    igl::gaussian_curvature(V,F,K);
    // Compute mass matrix
    Eigen::SparseMatrix<double> Minv;
    igl::invert_diag(M,Minv);
    // Divide by area to get integral average
    K = (Minv*K).eval();
    return K;
}

void Lap::harmonics(int k, Eigen::MatrixXd &eigVec, Eigen::VectorXd &eigVal)
{
    Eigen::MatrixXd eigVec_t;
    Eigen::VectorXd eigVal_t;
    if(!igl::eigs((-L).eval(),M,k,igl::EIGS_TYPE_SM,eigVec_t,eigVal_t))
    {
      std::cout<<"Lap::harmonics failed."<<std::endl;
    }
    // Normalize
    //eigVec_t = ((eigVec_t.array()-eigVec_t.minCoeff())/(eigVec_t.maxCoeff()-eigVec_t.minCoeff())).eval();
    eigVec.resize(eigVec_t.rows(), eigVec_t.cols());
    eigVal.resize(eigVal_t.size());
    for (int i = 0 ; i < k ; i++)
    {
        eigVec.col(i) = eigVec_t.col(k-1-i);
        eigVal(i) = eigVal_t(k-1-i);
    }
}

Eigen::SparseMatrix<double> Lap::getL() {return L;}
Eigen::SparseMatrix<double> Lap::getM() {return M;}
