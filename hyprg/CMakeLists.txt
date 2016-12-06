cmake_minimum_required(VERSION 3.6)
project(HypRG)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

set(SOURCE_FILES lib/graph.h lib/utils.h lib/graph.cpp lib/utils.cpp main.cpp lib/pair_generator.cpp lib/pair_generator.h lib/embedding_model.cpp lib/embedding_model.h lib/hyperbolic.cpp lib/hyperbolic.h lib/loss_function.cpp lib/loss_function.h lib/optimization.cpp lib/optimization.h)
add_executable(hyprg ${SOURCE_FILES})