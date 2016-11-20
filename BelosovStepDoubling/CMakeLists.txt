cmake_minimum_required(VERSION 3.6)
project(BelosovStepDoubling)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp stiff.cpp matrix.cpp utility.cpp)
add_executable(BelosovStepDoubling ${SOURCE_FILES})