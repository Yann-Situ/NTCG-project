/* laplacian.h
 * Author : Yann-Situ GAZULL
 * Description :
 */
#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

class Lap
{
private:
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    int n;
    Eigen::SparseMatrix<double> M; // mass matrix
    Eigen::SparseMatrix<double> L; // cotan Lap matrix
public:
    Lap(Eigen::MatrixXd _V, Eigen::MatrixXi _F);
    ~Lap();
    // Eigen::VectorXd curvature();
    void harmonics(int k, Eigen::MatrixXd &eigVec, Eigen::VectorXd &eigVal);
    void harmonics2(int k, Eigen::MatrixXd &eigVec, Eigen::VectorXd &eigVal);

    Eigen::SparseMatrix<double> getL();
    Eigen::SparseMatrix<double> getM();
};


#endif
