cmake_minimum_required(VERSION 3.6)
project(HypRG)

set (CMAKE_CXX_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lboost_program_options")
link_directories(/usr/local/lib)
include_directories(/usr/local/lib)

set(SOURCE_FILES main.cpp lib/graph.h lib/utils.h lib/graph.cpp lib/utils.cpp)
add_executable(hyprg ${SOURCE_FILES})
target_link_libraries(hyprg boost_program_options)