/* MHT.cpp
 * Author : Yann-Situ GAZULL
 * Description :
 */
#include "MHT.h"

#include <assert.h>
#include <math.h>       /* sqrt */
#include <iostream>


HarmonicTransform::HarmonicTransform(){}
HarmonicTransform::~HarmonicTransform(){}

Eigen::MatrixXd HarmonicTransform::MHT(Eigen::MatrixXd V,
    Eigen::SparseMatrix<double> M,
    Eigen::MatrixXd eigVec)
{
    int n = V.rows();  // nb vertices
    int m = eigVec.cols(); // nb harmonics
    assert(n == M.rows() && n == M.cols());
    assert(n == eigVec.rows());
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(m, V.cols());
    for (int d = 0 ; d < V.cols() ; d++)
    {
        R.col(d) = (V.col(d).transpose() * M * eigVec).transpose();
    }
    std::cout << "R :\n" << R.block(0,0,std::min(m,20),3) << "\n";
    return R;
}

Eigen::MatrixXd HarmonicTransform::invMHT(Eigen::MatrixXd R,
    Eigen::SparseMatrix<double> M,
    Eigen::MatrixXd eigVec)
{
    int m = R.rows();  // nb vertices
    int n = eigVec.rows(); // nb harmonics
    assert(n == M.rows() && n == M.cols());
    assert(m == eigVec.cols());
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, R.cols());
    for (int d = 0 ; d < R.cols() ; d++)
    {
        V.col(d) = eigVec * R.col(d);
        // for (int k = 0 ; k < m ; k++)
        // {
        //     V.col(d) = V.col(d) + R(k,d) * eigVec.col(k);
        // }

        // for (int i = 0 ; i < n ; i++)
        // {
        //     V(i,d) *= sqrt(M.coeff(i,i));
        // }
    }
    return V;
}
