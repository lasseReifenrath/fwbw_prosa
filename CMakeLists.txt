cmake_minimum_required(VERSION 3.27)
project(fwbw)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

set(CMAKE_CXX_STANDARD 17)
include(FindASan)
include_directories(simde)
include_directories(cacode)
add_subdirectory(cacode)

add_library(MyLib $<TARGET_OBJECTS:cacode> BaseMatrix.cpp SubstitutionMatrix.cpp FwBwAligner.cpp)
add_executable(fwbw main.cpp)
target_link_libraries(fwbw PRIVATE MyLib c++abi c++)