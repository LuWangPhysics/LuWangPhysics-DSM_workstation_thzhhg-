#description     :install dsm and required external libraries (fftw, eigen)
#author		 :Lu Wang adapted from Yikai Kan
#date            :2020/12/24
#==============================================================================

cmake_minimum_required(VERSION 3.0)
project(DSM_THz_hhg_LuWang)
include(ExternalProject)

add_custom_target(thirdparty)

ExternalProject_Add(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
    PREFIX ${PROJECT_SOURCE_DIR}/build/thirdparty/eigen
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${PROJECT_SOURCE_DIR}/build/thirdparty/eigen/ -DBUILD_SHARED_LIBS=OFF
)

ExternalProject_Add(
    fftw
    URL http://www.fftw.org/fftw-3.3.8.tar.gz
    PREFIX ${PROJECT_SOURCE_DIR}/build/thirdparty/fftw
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${PROJECT_SOURCE_DIR}/build/thirdparty/fftw/ -DBUILD_SHARED_LIBS=OFF
)

#ExternalProject_Add(
#    python
#    URL https://www.python.org/ftp/python/3.8.0/Python-3.8.0.tgz
#    PREFIX ${PROJECT_SOURCE_DIR}/build/thirdparty/python
#    CONFIGURE_COMMAND ./configure 
#	--enable-optimizations 
#	--prefix=${PROJECT_SOURCE_DIR}/build/thirdparty/python/  
#   BUILD_COMMAND make
#    BUILD_IN_SOURCE 1
#    INSTALL_COMMAND make install
#)



add_dependencies(thirdparty eigen fftw)


include_directories(${PROJECT_SOURCE_DIR}/build/thirdparty/eigen/)
include_directories(${PROJECT_SOURCE_DIR}/build/thirdparty/fftw/)
#include_directories(${PROJECT_SOURCE_DIR}/build/thirdparty/python/)

