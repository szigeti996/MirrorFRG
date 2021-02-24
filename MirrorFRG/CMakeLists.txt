# CMakeList.txt : CMake project for MirrorFRG, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)
project(Unitsv1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)  
# Add source to this project's executable.
set(SOURCE_FILES MirrorFRG.cpp) 
add_executable(Unitsv1 ${SOURCE_FILES})

# TODO: Add tests and install targets if needed.

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(Unitsv1 GSL::gsl GSL::gslcblas)