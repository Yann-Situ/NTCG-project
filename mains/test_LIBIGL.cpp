#include <igl/eigs.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/parula.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <iostream>
#include <queue>
#include "laplacian.h"

Eigen::MatrixXd V,U;
Eigen::MatrixXi F;
int c=0;
double bbd = 1;
bool twod = 0;
int main(int argc, char * argv[])
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  VectorXd D;
  const char* filename = (argc > 1) ? argv[1] : "../ressources/bunny.off";
  igl::readOFF(filename, V, F);


  Lap lap(V,F);
  twod = V.col(2).minCoeff()==V.col(2).maxCoeff();
  bbd = (V.colwise().maxCoeff()-V.colwise().minCoeff()).norm();
  SparseMatrix<double> L,M;
  cotmatrix(V,F,L);
  L = (-L).eval();
  massmatrix(V,F,MASSMATRIX_TYPE_DEFAULT,M);
  const size_t k = 5;

  lap.harmonics(k,U,D);

  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = [&](igl::opengl::glfw::Viewer & viewer,unsigned char key,int)->bool
  {
    switch(key)
    {
      default:
        return false;
      case ' ':
      {
        U = U.rightCols(k).eval();
        // Rescale eigen vectors for visualization
        VectorXd Z =
          bbd*0.5*U.col(c);
        if(twod)
        {
          V.col(2) = Z;
          viewer.data().set_mesh(V,F);
          viewer.data().compute_normals();
        }
        viewer.data().set_data(U.col(c).eval());
        c = (c+1)%U.cols();
        return true;
      }
    }
  };
  viewer.data().set_mesh(V,F);
  viewer.callback_key_down(viewer,' ',0);
  viewer.data().show_lines = false;
  std::cout<<
  R"(
  [space] Cycle through eigen modes
  )";
  std::cout<<D<<endl;;
  viewer.launch();
}

// #include <igl/gaussian_curvature.h>
// #include <igl/massmatrix.h>
// #include <igl/invert_diag.h>
// #include <igl/readOFF.h>
// #include <igl/opengl/glfw/Viewer.h>
// #include <iostream>
//
// #include "laplacian.h"
//
// int main(int argc, char *argv[])
// {
//   using namespace Eigen;
//   using namespace std;
//   MatrixXd V;
//   MatrixXi F;
//   const char* filename = (argc > 1) ? argv[1] : "../ressources/bunny.off";
//   igl::readOFF(filename, V, F);
//
//   Lap lap;
//   VectorXd K = lap.cuvature(V, F);
//   std::cout << K << "\n";
//
//   // Plot the mesh with pseudocolors
//   igl::opengl::glfw::Viewer viewer;
//   viewer.data().set_mesh(V, F);
//   viewer.data().set_data(K, -0.3, 0.3);
//   viewer.launch();
// }






// #include <igl/opengl/glfw/Viewer.h> // needs a link in the CMAKE
//
// int main(int argc, char *argv[])
// {
//   // Inline mesh of a cube
//   const Eigen::MatrixXd V= (Eigen::MatrixXd(8,3)<<
//     0.0,0.0,0.0,
//     0.0,0.0,1.0,
//     0.0,1.0,0.0,
//     0.0,1.0,1.0,
//     1.0,0.0,0.0,
//     1.0,0.0,1.0,
//     1.0,1.0,0.0,
//     1.0,1.0,1.0).finished();
//   const Eigen::MatrixXi F = (Eigen::MatrixXi(12,3)<<
//     1,7,5,
//     1,3,7,
//     1,4,3,
//     1,2,4,
//     3,8,7,
//     3,4,8,
//     5,7,8,
//     5,8,6,
//     1,5,6,
//     1,6,2,
//     2,6,8,
//     2,8,4).finished().array()-1;
//
//   // Plot the mesh
//   igl::opengl::glfw::Viewer viewer;
//   viewer.data().set_mesh(V, F);
//   viewer.data().set_face_based(true);
//   viewer.launch();
// }
