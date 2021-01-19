/* laplacian.cpp
 * Author : Yann-Situ GAZULL
 * Description :
 */
#include "laplacian.h"

//#include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/eigs.h>

//#include <Eigen/SparseCholesky>

#include <iostream>
#include <limits>

Lap::Lap(Eigen::MatrixXd _V, Eigen::MatrixXi _F)
: V(_V), F(_F)
{
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
    igl::cotmatrix(V,F,L);
    n = V.rows();
}
Lap::~Lap(){}

// Eigen::VectorXd Lap::curvature()
// {
//     Eigen::VectorXd K;
//     // Compute integral of Gaussian curvature
//     igl::gaussian_curvature(V,F,K);
//     // Compute mass matrix
//     Eigen::SparseMatrix<double> Minv;
//     igl::invert_diag(M,Minv);
//     // Divide by area to get integral average
//     K = (Minv*K).eval();
//     return K;
// }

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
    eigVec.resize(n, eigVec_t.cols());
    eigVal.resize(eigVal_t.size());
    for (int i = 0 ; i < k ; i++)
    {
        eigVec.col(i) = eigVec_t.col(k-1-i);
        eigVal(i) = eigVal_t(k-1-i);
    }
}

void Lap::harmonics2(int k, Eigen::MatrixXd &eigVec, Eigen::VectorXd &eigVal)
{
    eigVec.resize(n, k);
    eigVal.resize(k);
    double ls = 0;
    double llast = 0;
    double lk = 0;
    int i = 0; int k_temp = 0;
    Eigen::SparseMatrix<double> Minv;
    igl::invert_diag(M,Minv);
    Minv.cwiseSqrt();
    //Eigen::SparseMatrix<double> Lap = Minv * (-L).eval() * Minv;
    Eigen::SparseMatrix<double> Lap = (-L).eval();
    Eigen::SparseMatrix<double> Ds(n,n);
    Eigen::MatrixXd eigVec_t;
    Eigen::VectorXd eigVal_t;
    Eigen::SparseMatrix<double> I(n,n);
    I.setIdentity();
    std::cout << "GO FOR THE WHILE\n";

    while (i < k) {
        std::cout << " | NEW ITERATION ! : ";
        std::cout << i << "\n";
        Ds.setIdentity();
        Ds = Lap - ls * M;
        // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        // std::cout << i << "\n";
        // solver.compute(Ds);
        // std::cout << i << "\n";
        // Ds = solver.solve(I); // Ds^-1
        // std::cout << i << "\n";

        k_temp = std::min(k-i, 10);
        if(!igl::eigs(Ds,M, k_temp ,igl::EIGS_TYPE_SM,eigVec_t,eigVal_t))
        {
          std::cout << "Lap::harmonics " << i << " failed." << std::endl;
        }
        double maxi = 0.0;
        double mini = std::numeric_limits<double>::max();
        for (int j = 0; j < k_temp; j++) {
            // lk = ls+1.0/(eigVal_t(k_temp-1-j));
            lk = ls+eigVal_t(k_temp-1-j);
            if (lk > llast)
            {
                eigVec.col(i+j) = eigVec_t.col(k_temp-1-j);
                eigVal(i+j) = lk;
                maxi = std::max(maxi, lk);
                mini = std::min(mini, lk);
            }
        }
        i+=k_temp;
        std::cout << i << "\n";
        ls = maxi + 0.4 * (maxi - mini);
        llast = maxi;
    }
    //eigVec = Minv * eigVec;
}

Eigen::SparseMatrix<double> Lap::getL() {return L;}
Eigen::SparseMatrix<double> Lap::getM() {return M;}
