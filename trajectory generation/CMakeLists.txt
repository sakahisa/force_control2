cmake_minimum_required(VERSION 2.8)
project (HELLO)
#set(CMAKE_CXX_FLAGS "-Weffc++ -Wuninitialized")
#set(CMAKE_CXX_FLAGS "-fpermissive")

set(_srcs
	t-gene.cpp	)


#include_directories (Eigen)
include_directories("${PROJECT_SOURCE_DIR}")

add_executable (helloDemo ${_srcs})

set(CMAKE_BUILD_TYPE Debug)
