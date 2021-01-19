/* filter.h
 * Author : Yann-Situ GAZULL
 * Description :
 */
#ifndef FILTER_H
#define FILTER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>


class Filter
{
private:
    Filter();
    ~Filter();


public:

    static Eigen::MatrixXd apply(Eigen::MatrixXd V, Eigen::MatrixXd R,
        Eigen::MatrixXd eigVec,
        Eigen::SparseMatrix<double> M,
        Eigen::VectorXd F,
        double filter_limit_value = 1.0);

    static Eigen::MatrixXd applyNoise(Eigen::MatrixXd V, Eigen::MatrixXi F,
         double amp = 0.05);

    /*
    * Filtering high frequancies with F(x) = 1 if x < w_m Random
    * F(x) = lim + (1-lim)/(1+s_m(x-w_m)^2) otherwise
    */
    static Eigen::VectorXd inverseFilter(double w_m, double s_m, double lim,
        Eigen::VectorXd freq);

    /*
    * Filtering high frequencies with
    * F(x) = 1 + (value-1)/(1+s_m(x-w_m)^2) otherwise
    */
    static Eigen::VectorXd focusedFilter(double w_m, double s_m, double value,
        Eigen::VectorXd freq);

    static Eigen::VectorXd constantFilter(double val, int freq_size);
};



#endif
