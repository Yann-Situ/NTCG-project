/* MHT.h
 * Author : Yann-Situ GAZULL
 * Description : Manifold Harmonic transform
 */
#ifndef MHT_H
#define MHT_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

class HarmonicTransform
{
private:


public:
    HarmonicTransform();
    ~HarmonicTransform();
    /*
    * Perform a Manifold Harmonic transform, in order to get the coeff of
    * harmonic decomposition from :
    * V   : the vertex positions
    * M   : the mass matrix / Hodge star
    * eigVec : the eigen vectors of the symetrized laplacian
    *          (M^{-1/2}LM^{1/2} where L is the Beltrami with cotans)
    * return a vector of the eigVec.cols() coeff of the harmonic decomposition.
    */
    Eigen::MatrixXd MHT(Eigen::MatrixXd V,
        Eigen::SparseMatrix<double> M,
        Eigen::MatrixXd eigVec);


    /*
    * Perform the inverse Manifold Harmonic transform, in order to get the
    * positions of the vertices from the harmonic decomposition.
    * R   : the harmonic decomposition coefficients
    * M   : the mass matrix / Hodge star
    * eigVec : the eigen vectors of the symetrized laplacian
    *          (M^{-1/2}LM^{1/2} where L is the Beltrami with cotans)
    * return the eigVec.rows() vertex positions
    */
    Eigen::MatrixXd invMHT(Eigen::MatrixXd R,
        Eigen::SparseMatrix<double> M,
        Eigen::MatrixXd eigVec);
};

#endif
