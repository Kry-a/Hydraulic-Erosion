cmake_minimum_required(VERSION 3.0)
project(Hydraulic-Erosion)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(include/random)

find_package(TIFF)
find_package(OpenMP)

include_directories("include/SimplexNoise/src")
add_library(SimplexNoise STATIC include/SimplexNoise/src/SimplexNoise.cpp)

add_executable(${PROJECT_NAME} src/main.cpp src/Erosion.cpp)
target_link_libraries(${PROJECT_NAME} effolkronium_random SimplexNoise TIFF::TIFF OpenMP::OpenMP_CXX)
