cmake_minimum_required(VERSION 3.1)

project(NTCG-project)

#set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#set(CMAKE_CXX_FLAGS "-O2") #slow a bit the compilation (source : https://www.cgal.org/FAQ.html#compilation_speed)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(COMPILE_TEST true)

###########################################

file(GLOB INCLUDES include/*.h)
file(GLOB SOURCES src/*.cpp)
add_library( BaseFiles ${SOURCES} ${INCLUDES} )

include_directories(
    include/
    include/libigl/include
)

###########################################

# Eigen
find_package(Eigen3 3.1.0 REQUIRED)
if ( EIGEN3_FOUND )
    include_directories(${EIGEN3_INCLUDE_DIR})
    #include( ${EIGEN3_USE_FILE} ) # is equivalent (maybe ?)
    set(CGAL_EIGEN3_ENABLED true)
    message(STATUS "Eigen3 found : " ${EIGEN3_INCLUDE_DIR})
else()
    message(SEND_ERROR "ERROR: this program requires the EIGEN3 library")
    return()
endif()

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON) # for the viewer only
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON) # for the viewer only

find_package(LIBIGL REQUIRED QUIET)

###########################################

# Create the executable of the application
# No need to link libraries for eigen3 and libigl (except for viewer)
# as they are header only.

if ( COMPILE_TEST )
    #EIGEN3 TEST
    set(exe-test_EIGEN3 "test_EIGEN3")
    add_executable(exe-test_EIGEN3  mains/test_EIGEN3.cpp)
    target_link_libraries(
        exe-test_EIGEN3
        BaseFiles
    )

    #LIBIGL TEST
    set(exe-test_LIBIGL "test_LIBIGL")
    add_executable(exe-test_LIBIGL  mains/test_LIBIGL.cpp)
    target_link_libraries(
        exe-test_LIBIGL
        BaseFiles
        igl::core igl::opengl_glfw
    )
endif()
