cmake_minimum_required(VERSION 3.10)
project(pyflagsercount VERSION 0.0.1)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# flagser-count-maximal
add_executable(flagser-count-maximal src/test_maximal.cpp)

#pthread
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(flagser-count-maximal PRIVATE Threads::Threads)


#flagser-count-filtered
include_directories(sparsehash/src)

add_executable(flagser-count-filtered src/test_filtered.cpp)
target_link_libraries(flagser-count-filtered PRIVATE Threads::Threads)


#flagser-contain binding with pybind11
add_subdirectory(pybind11)
pybind11_add_module(pyflagsercontain src/flagser_count_bindings.cpp)



