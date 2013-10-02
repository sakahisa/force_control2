cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(fastgetsetdcm)
find_package(qibuild)

set(_srcs
    main.cpp
    fast.h
    fast.cpp
)

qi_create_lib(fastgetsetdcm SHARED ${_srcs} SUBFOLDER naoqi)

qi_use_lib(fastgetsetdcm ALCOMMON ALMEMORYFASTACCESS)
