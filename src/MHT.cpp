/* MHT.cpp
 * Author : Yann-Situ GAZULL
 * Description :
 */
#include "MHT.h"

#include <assert.h>
#include <math.h>       /* sqrt */
#include <iostream>

Eigen::MatrixXd HarmonicTransform::MHT(Eigen::MatrixXd V,
    Eigen::MatrixXd M,
    Eigen::MatrixXd eigVec)
{
    int n = V.rows();  // nb vertices
    int m = eigVec.cols(); // nb harmonics
    assert(n == M.rows() && n == M.cols());
    assert(n == eigVec.rows());

    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(m, V.cols());
    for (int d = 0 ; d < V.cols() ; d++)
    {
        for (int i = 0 ; i < n ; i++)
        {
            R.col(d) = R.col(d) + V(i,d) * sqrt(M(i,i)) * eigVec.row(i).tranpose();
        }
    }
    return R;
}

Eigen::MatrixXd HarmonicTransform::invMHT(Eigen::MatrixXd R,
    Eigen::MatrixXd M,
    Eigen::MatrixXd eigVec)
{
    int m = R.rows();  // nb vertices
    int n = eigVec.rows(); // nb harmonics
    assert(n == M.rows() && n == M.cols());
    assert(m == eigVec.cols());

    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, R.cols());
    for (int d = 0 ; d < R.cols() ; d++)
    {
        for (int i = 0 ; i < n ; i++)
        {
            for (int k = 0 ; k < m ; k++)
            {
                V(i,d) = V(i,d) + R(k) * sqrt(M(i,i)) * eigVec(i,k);
            }
        }
    }
    return V;
}
