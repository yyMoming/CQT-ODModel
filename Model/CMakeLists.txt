cmake_minimum_required(VERSION 3.18)
set(CMAKE_CXX_STANDARD 14)
project(Model)
#set(Torch_DIR /home/data/ywm_data/libtorch/share/cmake/Torch)
#find_package(Torch REQUIRED)
#if(TORCH_FOUND)
#    message(status "@@@@@@@@@@ torch has been included")
#else()
#    message(status "@@@@@@@@@@ torch not found")
#endif(TORCH_FOUND)
#
#message(status "@@@@@@@@@@ " ${PROJECT_SOURCE_DIR})
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

file(GLOB HXX include/*)
aux_source_directory(src CXX)   #只收集源文件


message("@@@@@@@@@@ headers: " ${HXX} )
message("@@@@@@@@@@ sources: " ${CXX} )
add_library(${PROJECT_NAME} ${HXX} ${CXX})
#add_executable(${PROJECT_NAME} ${HXX} ${CXX})
target_include_directories(${PROJECT_NAME} PUBLIC ./include)
#message("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 1")
target_link_libraries(${PROJECT_NAME} PUBLIC "${TORCH_LIBRARIES}")
#message("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 2")
#target_include_directories(${PROJECT_NAME} PUBLIC ../utility)

