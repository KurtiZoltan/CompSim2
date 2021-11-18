cmake_minimum_required(VERSION 3.5)

project(CompSim1)

set(CMAKE_CXX_STANDARD 17) #for aligned_alloc
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_VERBOSE_MAKEFILE OFF)

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    add_compile_options(/W4)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-Wall -Wextra -pedantic -Wfatal-errors -fvisibility=hidden)
endif()

file(GLOB SOURCE "src/*.cpp" "src/*.hpp")

find_package(TIFF REQUIRED)
#find_package(OpenMP REQUIRED)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

add_executable(1d64-32laplace ./src/1d64-32laplace.cpp)
add_executable(2dlaplacetest ./src/2dlaplacetest.cpp)
add_executable(timeEvolution ./src/timeEvolution.cpp)
