cmake_minimum_required(VERSION 3.11)


project(x2chem)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES x2chem.cpp x2chem.hpp)


# CMake instructions to make the static lib
add_library(x2chem STATIC ${SOURCE_FILES})



# CMake instructions to test using the static lib
set( APP_EXE StaticTest )

add_executable( ${APP_EXE} main.cxx ) 

target_link_libraries( ${APP_EXE} x2chem )
