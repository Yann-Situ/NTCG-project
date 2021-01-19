# NTCG-project

- *LIBIGL* (https://github.com/libigl/libigl.git) must be cloned in `./include`.
- If it is already on your computer, it can be linked to the project by modifying the `CMakLists.txt` using `include_directories(PATH_TO_LIBIGL/include)`.

### Dependencies

The only dependencies are stl, eigen, [libigl](http://libigl.github.io/libigl/) and
the dependencies of the `igl::opengl::glfw::Viewer`.

## Manifold Harmonics
### Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

### Run

Run `exe-main` after compilation.
A glfw app should launch *ressources/bunny.off*.
You can launch your own mesh with `exe-main my_own_mesh.off`.
