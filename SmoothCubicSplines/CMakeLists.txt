cmake_minimum_required(VERSION 3.6)
project(SmoothCubicSplines)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp hermiteCubic.cpp hermiteCubic.h matrix.cpp matrix.h utility.cpp utility.h)
add_executable(SmoothCubicSplines ${SOURCE_FILES})