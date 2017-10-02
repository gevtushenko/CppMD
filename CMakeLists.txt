cmake_minimum_required(VERSION 3.8)
project(MD)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

set(SOURCE_FILES main.cpp)
add_executable(MD ${SOURCE_FILES})
