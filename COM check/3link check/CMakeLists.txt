cmake_minimum_required(VERSION 2.8)
project (HELLO)
#set(CMAKE_CXX_FLAGS "-Weffc++ -Wuninitialized")
set(CMAKE_CXX_FLAGS "-fpermissive")

#include_directories (Eigen)
include_directories("${PROJECT_SOURCE_DIR}")

add_executable (helloDemo main.cpp)

set(CMAKE_BUILD_TYPE Debug)
