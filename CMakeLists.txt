cmake_minimum_required(VERSION 3.12)
project(MathOMP VERSION 0.1.0)
include(GNUInstallDirs)
SET(CMAKE_CXX_FLAGS "-lstdc++ -O3 -Wall -pthread -fcheck-new -g  ")
set(CMAKE_CXX_STANDARD 14)
find_package(OpenMP REQUIRED)

link_libraries(-lm -lpthread)

include_directories("include/mathomp" "include" "/Users/navgup/repos/benchmark/include")

SET(GOOGLE_BENCHMARK_SRC src/gbenchmark.cpp)



set(SOURCE_FILES src/mathomp.cpp include/mathomp/mathomp.h)




add_library(MathOMP STATIC ${SOURCE_FILES})


if (APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "$(brew --prefix libomp)/lib/libomp.dylib")
endif()

if(OpenMP_CXX_FOUND)
endif()

install(TARGETS MathOMP DESTINATION lib)
install(FILES mathomp.h DESTINATION include)
