cmake_minimum_required(VERSION 3.6)
project(GalerkinMethod)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp gal2ptBVP.cpp)
add_executable(GalerkinMethod ${SOURCE_FILES})