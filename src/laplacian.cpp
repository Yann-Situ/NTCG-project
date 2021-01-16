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
    if(!igl::eigs(L,M,k,igl::EIGS_TYPE_SM,eigVec,eigVal))
    {
      std::cout<<"Lap::harmonics failed."<<std::endl;
    }
    // Normalize
    eigVec = ((eigVec.array()-eigVec.minCoeff())/(eigVec.maxCoeff()-eigVec.minCoeff())).eval();
}

Eigen::SparseMatrix<double> Lap::getL() {return L;}
Eigen::SparseMatrix<double> Lap::getM() {return M;}
