cmake_minimum_required(VERSION 3.18)
set(CMAKE_CXX_STANDARD 14)
project(utility)
aux_source_directory(. CXX)
add_library(${PROJECT_NAME} ${CXX})