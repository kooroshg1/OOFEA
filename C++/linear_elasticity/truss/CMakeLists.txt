cmake_minimum_required(VERSION 3.2)
project(vulcan-development)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# ======================================================================== #
# INCLUDE Directories
# find libmesh include
include_directories($ENV{LIBMESH_DIR}/include)
# find mpi include
include_directories(/usr/local/include)
# find extra includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
# ======================================================================== #
file(GLOB SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
        )
# ======================================================================== #
# LIB Directories
# find libmesh library
find_library(libmesh-lib
        NAMES mesh_opt
        HINTS $ENV{LIBMESH_DIR}/lib
        NO_DEFAULT_PATH NO_CMAKE_FIND_ROOT_PATH)
# find mpi library
find_library(mpi-lib
        NAMES mpi
        HINTS /usr/local/lib
        NO_DEFAULT_PATH NO_CMAKE_FIND_ROOT_PATH)
find_library (mpicxx_lib
        NAMES mpi_cxx
        HINTS /usr/local/lib
        NO_DEFAULT_PATH NO_CMAKE_FIND_ROOT_PATH)
# ======================================================================== #
find_package (Threads REQUIRED)
set(SOURCE_FILES main.cpp)
add_executable(vulcan ${SOURCE_FILES} ${SOURCES})
TARGET_LINK_LIBRARIES(vulcan
        ${libmesh-lib}
        ${mpi-lib}
        ${mpicxx_lib}
        -lpthread ${CMAKE_THREAD_LIBS_INIT})