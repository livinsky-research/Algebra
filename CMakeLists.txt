cmake_minimum_required(VERSION 3.5)

project(Algebra)

set(CMAKE_VERBOSE_MAKEFILE off)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#set(CMAKE_CXX_FLAGS "-lstdc++")
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pthread")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
#set(CMAKE_BUILD_TYPE ASAN)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

find_package(Catch REQUIRED)

set(SOURCES
    src/Permutation.cpp
    src/Group.cpp
    src/Structure.cpp
)

include_directories(
    ${PROJECT_SOURCE_DIR}/inc
    ${PROJECT_SOURCE_DIR}/catch
)

add_subdirectory(tests)
add_subdirectory(examples)
