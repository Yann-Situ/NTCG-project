#NTCG-project

- *LIBIGL* (https://github.com/libigl/libigl.git) must be cloned in `./include`.
- If it is already on your computer, it can be linked to the project by modifying the `CMakLists.txt` using `include_directories(PATH_TO_LIBIGL/include)`.

## libigl example project

A blank project example showing how to use libigl and cmake. Feel free and
encouraged to copy or fork this project as a way of starting a new personal
project using libigl.

### See the tutorial first

Then build, run and understand the [libigl
tutorial](http://libigl.github.io/libigl/tutorial/).

### Dependencies

The only dependencies are stl, eigen, [libigl](http://libigl.github.io/libigl/) and
the dependencies of the `igl::opengl::glfw::Viewer`.

The cmake build system will attempt to find libigl according to environment variables (e.g., `LIBIGL`) and searching in common desitinations (e.g., `/usr/local/libigl/`). If you haven't installed libigl before, we recommend you to clone a copy of libigl right here:

    cd libigl-example-project/
    git clone https://github.com/libigl/libigl.git

### Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `example_bin` binary.

### Run

From within the `build` directory just issue:

    ./example

A glfw app should launch displaying a 3D cube.

### Time
10 minutes on paper + 5 minute on project
Report : 1/3 on paper, 2/3 on project
