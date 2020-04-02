cmake_minimum_required(VERSION 3.15)
project(cpp)

set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
add_executable(cpp main.cpp integrand_params.cpp integrand_params.h)
find_package(GSL REQUIRED)
#find_package(OpenMP REQUIRED)
target_link_libraries(cpp GSL::gsl GSL::gslcblas)