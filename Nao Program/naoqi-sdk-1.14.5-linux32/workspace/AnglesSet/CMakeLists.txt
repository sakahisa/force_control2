cmake_minimum_required(VERSION 2.8)
project(AnglesSet)

find_package(qibuild)

#Create a list of source files
set(_srcs
        Set.h
        Set.cpp
	AnglesSetModule.cpp
	AnglesSetModule.h
	main.cpp)

include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/Eigen")

#Create a plugin, that is a shared library, and make
#sure it is built in lib/naoqi
qi_create_lib(AnglesSetModule SHARED ${_srcs} SUBFOLDER naoqi)

#Tell CMake that AnglesSetModule depends on ALCOMMON.
#This will set the libraries to link AnglesSetModule with,
#the include paths, and so on
qi_use_lib(AnglesSetModule ALCOMMON BOOST)
