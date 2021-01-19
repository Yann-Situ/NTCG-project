/* filter.cpp
 * Author : Yann-Situ GAZULL
 * Description :
 */
#include "MHT.h"
#include "filter.h"
#include <random>
#include <ctime>
#include <igl/per_vertex_normals.h>

Filter::Filter(){}
Filter::~Filter(){}

Eigen::MatrixXd Filter::apply(Eigen::MatrixXd V, Eigen::MatrixXd R,
    Eigen::MatrixXd eigVec,
    Eigen::SparseMatrix<double> M,
    Eigen::VectorXd F,
    double filter_limit_value)
{
    Eigen::MatrixXd Vp(V.rows(), V.cols());
    Eigen::MatrixXd Vhf=V-HarmonicTransform::invMHT(R, M, eigVec);
    Vp = HarmonicTransform::invMHT(F.asDiagonal() * R, M, eigVec);
    Vp += filter_limit_value * Vhf;
    return Vp;
}

Eigen::MatrixXd Filter::applyNoise(Eigen::MatrixXd V, Eigen::MatrixXi F,
     double amp)
{
    Eigen::MatrixXd Vp = V;
    Eigen::MatrixXd N;
    igl::per_vertex_normals(V,F,N);

    std::srand(std::time(nullptr));
    for (int i = 0; i < V.rows(); i++) {
        double r = -1.0 + 2.0*(double)(rand()) / ((double)(RAND_MAX));
        Vp.row(i) += r*amp*N.row(i);
    }
    return Vp;
}

Eigen::VectorXd Filter::inverseFilter(double w_m, double s_m, double lim,
    Eigen::VectorXd freq)
{
    Eigen::VectorXd F(freq.size());
    for (int i = 0; i < freq.size(); i++) {
        if (freq(i) > w_m)
        {
            double temp = (freq(i)-w_m);
            F(i) = lim + (1.0-lim)*1.0/(1.0+s_m*temp*temp);
        }
        else
        {
            F(i) = 1.0;
        }
    }
    return F;
}

Eigen::VectorXd Filter::focusedFilter(double w_m, double s_m, double value,
    Eigen::VectorXd freq)
{
    Eigen::VectorXd F(freq.size());
    for (int i = 0; i < freq.size(); i++) {
        double temp = (freq(i)-w_m);
        F(i) = 1.0 + (value-1.0)/(1.0+s_m*temp*temp);
    }
    return F;
}

Eigen::VectorXd Filter::constantFilter(double val, int freq_size)
{
    Eigen::VectorXd F(freq_size);
    for (int i = 0; i < freq_size; i++) {
        F(i) = val;
    }
    return F;
}
