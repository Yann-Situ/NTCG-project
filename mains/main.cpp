/* main.cpp
 * Author : Yann-Situ GAZULL
 * Description :
 */
#include <igl/readOFF.h>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/invert_diag.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

#include <igl/png/writePNG.h>

#include "laplacian.h"
#include "MHT.h"
#include "filter.h"

#include <iostream>
#include <sstream>
#include <string>

#define NB_HARMONICS_DEFAULT 30

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  Eigen::MatrixXd V; // vertices
  Eigen::MatrixXi F; // faces
  const char* filename = (argc > 1) ? argv[1] : "../ressources/bunny.off";
  igl::readOFF(filename, V, F);
  V = (-0.5+(V.array()-V.minCoeff())/(V.maxCoeff()-V.minCoeff())).eval();
  // Eigen::MatrixXd C_default_color(V.rows(),3);
  // for (int i = 0 ; i<C_default_color.rows() ; i++)
  // {
  //     C_default_color(i,0) = 0.9;
  //     C_default_color(i,1) = 0.8;
  //     C_default_color(i,2) = 0.3;
  // }

  Lap Lap_handler(V,F);
  //Eigen::VectorXd K = Lap_handler.curvature();
  //K = ((K.array()-K.minCoeff())/(K.maxCoeff()-K.minCoeff())).eval(); (normalize)
  //std::cout << K << "\n";
  Eigen::MatrixXd eigVec; Eigen::VectorXd eigVal; Eigen::VectorXd freq;
  Eigen::MatrixXd V_todisplay = V;
  Eigen::MatrixXd R;//harmonic coefficients

  // Customize the filtering menu
  double f_limit = 1.0;
  double f_flat = 0.0;
  double f_freq = 0.0;

  double noise_amp = 0.05;

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;

  // Attach a menu plugin
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);
  double doubleVariable = 1.0;
  int screenshot_count = 100;


  // Add content to the default menu window
  menu.callback_draw_viewer_menu = [&]()
  {
    // Draw parent menu content
    menu.draw_viewer_menu();
  };

  // Draw additional windows
  menu.callback_draw_custom_window = [&]()
  {
    // Define next window position + size
    ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(200, 400), ImGuiCond_FirstUseEver);
    ImGui::Begin(
        "Manifold Harmonics", nullptr,
        ImGuiWindowFlags_NoSavedSettings
    );

    // // Expose the same variable directly ...
    // ImGui::PushItemWidth(-80);
    // ImGui::DragScalar("double", ImGuiDataType_Double, &doubleVariable, 0.1, 0, 0, "%.2f");
    // ImGui::PopItemWidth();
    //
    // static std::string str = filename;
    // ImGui::InputText("Name", str);

    // Add new group
    // if (ImGui::CollapsingHeader("New Group", ImGuiTreeNodeFlags_DefaultOpen))
    // {
    //   // Expose variable directly ...
    //   ImGui::InputDouble("double", &doubleVariable, 0, 0, "%.4f");
    //
    //   // ... or using a custom callback
    //   static bool boolVariable = false;
    //   if (ImGui::Checkbox("bool", &boolVariable))
    //   {
    //     // do something
    //     std::cout << "boolVariable: " << std::boolalpha << boolVariable << std::endl;
    //   }
    //
    //   // Expose an enumeration type
    //   enum Orientation { Up=0, Down, Left, Right };
    //   static Orientation dir = Up;
    //   ImGui::Combo("Direction", (int *)(&dir), "Up\0Down\0Left\0Right\0\0");

      // // We can also use a std::vector<std::string> defined dynamically
      // static int num_choices = 3;
      // static std::vector<std::string> choices;
      // static int idx_choice = 0;
      // if (ImGui::InputInt("Num letters", &num_choices))
      // {
      //   num_choices = std::max(1, std::min(26, num_choices));
      // }
      // if (num_choices != (int) choices.size())
      // {
      //   choices.resize(num_choices);
      //   for (int i = 0; i < num_choices; ++i)
      //     choices[i] = std::string(1, 'A' + i);
      //   if (idx_choice >= num_choices)
      //     idx_choice = num_choices - 1;
      // }
      // ImGui::Combo("Letter", &idx_choice, choices);


      static int harmonics_nb = NB_HARMONICS_DEFAULT;
      static std::vector<std::string> harmonics_choices = {"no display"};
      static int idx_harmonics_choice = 0;
      if (ImGui::CollapsingHeader("Harmonics", ImGuiTreeNodeFlags_DefaultOpen))
      {
          if (ImGui::InputInt("Num harmonics to compute", &harmonics_nb))
          {
              harmonics_nb = std::max(1, std::min(1000, harmonics_nb));
          }
          // if (ImGui::Button("Compute harmonics", ImVec2(-1,0)))
          // {
          //     std::cout << "computing the " << harmonics_nb << " first harmonics...\n";
          //     Lap_handler.harmonics(harmonics_nb, eigVec, eigVal);
          //     freq = eigVal.cwiseSqrt();
          //     //Eigen::SparseMatrix<double> M = Lap_handler.getM();
          //     std::cout << "Done! Eigen values are " << eigVal << "\n";
          //     //std::cout << "eigVec :\n" << eigVec.block(0,0,10,5) << "\n";
          //     //std::cout << "eigVecT*M*eigVec :\n" << eigVec.transpose()* M * eigVec << "\n"; // should be identity
          // }
          if (ImGui::Button("Compute harmonics", ImVec2(-1,0)))
          {
              std::cout << "computing the " << harmonics_nb << " first harmonics...\n";
              Lap_handler.harmonics2(harmonics_nb, eigVec, eigVal);
              freq = eigVal.cwiseSqrt();
              f_freq = freq(harmonics_nb / 2);
              f_flat = 1.0/(freq(3*harmonics_nb / 4)-freq(harmonics_nb / 2));
              //Eigen::SparseMatrix<double> M = Lap_handler.getM();
              std::cout << "Done! Eigen values are " << eigVal << "\n";
              //std::cout << "eigVec :\n" << eigVec.block(0,0,10,5) << "\n";
              //std::cout << "eigVecT*M*eigVec :\n" << eigVec.transpose()* M * eigVec << "\n"; // should be identity
          }

          if (eigVec.cols()+1 != (int) harmonics_choices.size())
          {
              harmonics_choices.resize(eigVec.cols()+1);
              for (int i = 1; i < eigVec.cols()+1; ++i)
              harmonics_choices[i] = std::to_string(i);
              if (idx_harmonics_choice > eigVec.cols())
              idx_harmonics_choice = eigVec.cols();
          }
          //ImGui::Combo("Harmonic to display", &idx_harmonics_choice, harmonics_choices);
          if (ImGui::Combo("Display harmonic", &idx_harmonics_choice, harmonics_choices))
          {
              if (idx_harmonics_choice >= 1 && eigVec.cols() >= idx_harmonics_choice)
              {
                  viewer.data().set_data(eigVec.col(idx_harmonics_choice-1).eval());
                  std::cout << "Harmonic num " << idx_harmonics_choice << "\n";
                  std::cout << " | eigen value : " << eigVal(idx_harmonics_choice-1) << "\n";
                  std::cout << " | max value   : " << eigVec.col(idx_harmonics_choice-1).maxCoeff() << "\n";
                  std::cout << " | min value   : " << eigVec.col(idx_harmonics_choice-1).minCoeff() << "\n";
              }
              else
              {
                  std::cout << "Not displaying harmonics\n";
                  viewer.data().set_data(Filter::constantFilter(1.0,V.rows()));
                  // viewer.data().set_colors(C_default_color);
              }
          }
      }
      if (ImGui::CollapsingHeader("Filtering", ImGuiTreeNodeFlags_DefaultOpen))
      {

          if (ImGui::Button("F : forget HF", ImVec2(-1,0)))
          {
              R = HarmonicTransform::MHT(V, Lap_handler.getM(), eigVec);
              std::cout << "R dims : " << R.rows() << "x"
              << R.cols() << std::endl;
              V_todisplay = Filter::apply(V, R, eigVec, Lap_handler.getM(),
              Filter::constantFilter(1.0,freq.size()), 0.0);
              //V_todisplay = HarmonicTransform::invMHT(R, Lap_handler.getM(), eigVec);
              viewer.data().set_vertices(V_todisplay);
              std::cout << "Filter applied\n";
          }

          if (ImGui::Button("F : forget LF", ImVec2(-1,0)))
          {
              R = HarmonicTransform::MHT(V, Lap_handler.getM(), eigVec);
              std::cout << "R dims : " << R.rows() << "x"
              << R.cols() << std::endl;
              V_todisplay = Filter::apply(V, R, eigVec, Lap_handler.getM(),
              Filter::constantFilter(0.0,freq.size()), 1.0);
              //V_todisplay = HarmonicTransform::invMHT(R, Lap_handler.getM(), eigVec);
              viewer.data().set_vertices(V_todisplay);
              std::cout << "Filter applied\n";
          }

          ImGui::InputDouble("filter limit/value", &f_limit, 0, 0, "%.2f");
          ImGui::InputDouble("filter flatness", &f_flat, 0, 0, "%.2f");
          ImGui::InputDouble("filter thresh frequency", &f_freq, 0, 0, "%.2f");

          if (ImGui::Button("F : asymptotic filtering", ImVec2(-1,0)))
          {
              R = HarmonicTransform::MHT(V, Lap_handler.getM(), eigVec);
              std::cout << "R dims : " << R.rows() << "x"
                                       << R.cols() << std::endl;

              Eigen::VectorXd filter;
              filter = Filter::inverseFilter(f_freq, f_flat, f_limit, freq);
              if (f_limit > 1.0)
                  std::cout << " | filtering low frequencies\n";
              else
                  std::cout << " | filtering high frequencies\n";

              V_todisplay = Filter::apply(
                  V, R, eigVec, Lap_handler.getM(),
                  filter, f_limit);
              //V_todisplay = HarmonicTransform::invMHT(R, Lap_handler.getM(), eigVec);
              viewer.data().set_vertices(V_todisplay);
              std::cout << "Filter applied\n";
          }
          if (ImGui::Button("F : focused filtering", ImVec2(-1,0)))
          {
              R = HarmonicTransform::MHT(V, Lap_handler.getM(), eigVec);
              std::cout << "R dims : " << R.rows() << "x"
                                       << R.cols() << std::endl;

              Eigen::VectorXd filter;
              filter = Filter::focusedFilter(f_freq, f_flat, f_limit, freq);

              V_todisplay = Filter::apply(
                  V, R, eigVec, Lap_handler.getM(),
                  filter, 1.0);
              //V_todisplay = HarmonicTransform::invMHT(R, Lap_handler.getM(), eigVec);
              viewer.data().set_vertices(V_todisplay);
              std::cout << "Filter applied\n";
          }
      }

      if (ImGui::CollapsingHeader("Misc", ImGuiTreeNodeFlags_DefaultOpen))
      {
          if (ImGui::Button("Reset vertices", ImVec2(-1,0)))
          {
              V_todisplay = V;
              viewer.data().set_vertices(V_todisplay);
          }
          ImGui::InputDouble("noise amplitude", &noise_amp, 0, 0, "%.3f");
          if (ImGui::Button("Apply noise", ImVec2(-1,0)))
          {
              V_todisplay = Filter::applyNoise(V, F, noise_amp);
              viewer.data().set_vertices(V_todisplay);
          }
          if (ImGui::Button("ScreenShot", ImVec2(-1,0)))
          {
              viewer.data().set_mesh(V_todisplay, F);
              //viewer.data().set_data(K);
              //viewer.data().add_label(viewer.data().V.row(0) + viewer.data().V_normals.row(0).normalized()*0.005, "Hello World!");

              //viewer.launch();

              // Launch viewer and draw
              //viewer.launch_init();
              //viewer.draw();
              // Allocate temporary buffers
              Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Im_R(1280,800);
              Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Im_G(1280,800);
              Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Im_B(1280,800);
              Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Im_A(1280,800);

              // Draw the scene in the buffers
              viewer.core().draw_buffer(viewer.data(),false,Im_R,Im_G,Im_B,Im_A);

              // Save it to a PNG
              std::stringstream result_name;
              result_name << "../results/out_";
              result_name << screenshot_count << ".png";
              screenshot_count++;
              igl::png::writePNG(Im_R,Im_G,Im_B,Im_A, result_name.str());
              std::cout << "ScreenShot " << result_name.str() << "\n";

              //viewer.launch_shut();
          }
      }
      // if (ImGui::Button("F : keep all", ImVec2(-1,0)))
      // {
      //     R = HarmonicTransform::MHT(V, Lap_handler.getM(), eigVec);
      //     std::cout << "R dims : " << R.rows() << "x"
      //     << R.cols() << std::endl;
      //     V_todisplay = Filter::apply(V, R, eigVec, Lap_handler.getM(),
      //     Filter::constantFilter(1.0,freq), 1.0);
      //     //V_todisplay = HarmonicTransform::invMHT(R, Lap_handler.getM(), eigVec);
      //     viewer.data().set_vertices(V_todisplay);
      //     std::cout << "Filter applied\n";
      // }
    ImGui::End();
  };
  // Plot the mesh
  viewer.data().set_mesh(V_todisplay, F);
  //viewer.data().set_data(K);
  //viewer.data().add_label(viewer.data().V.row(0) + viewer.data().V_normals.row(0).normalized()*0.005, "Hello World!");

  viewer.launch();


}
