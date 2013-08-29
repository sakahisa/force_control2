project (HELLO)
#set(CMAKE_CXX_FLAGS "-Weffc++ -Wuninitialized")
set(CMAKE_CXX_FLAGS "-fpermissive")

#include_directories (Eigen)

add_executable (helloDemo main.cpp)

set(CMAKE_BUILD_TYPE Debug)


