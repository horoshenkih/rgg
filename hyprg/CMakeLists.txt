cmake_minimum_required(VERSION 3.6)
project(HypRG)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

set(SOURCE_FILES src/graph.h src/utils.h src/graph.cpp src/utils.cpp main.cpp src/pair_generator.cpp src/pair_generator.h src/embedding_model.cpp src/embedding_model.h src/hyperbolic.cpp src/hyperbolic.h src/loss_function.cpp src/loss_function.h src/optimization.cpp src/optimization.h src/commandline.h)
add_executable(hyprg ${SOURCE_FILES})