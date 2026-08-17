cmake_minimum_required(VERSION 3.16)
project(hp_fem CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

option(FEM_BUILD_TESTS "Build the test suite" ON)

find_package(Eigen3 REQUIRED NO_MODULE)

# everything except the driver lives in a library, so that both the executable
# and the tests can link against the same objects
file(GLOB FEM_SOURCES ${CMAKE_SOURCE_DIR}/source/*.cpp)

add_library(fem_core STATIC ${FEM_SOURCES})
target_include_directories(fem_core PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(fem_core PUBLIC Eigen3::Eigen)
target_compile_options(fem_core PRIVATE -Wall -Wextra)

add_executable(FEM ${CMAKE_SOURCE_DIR}/main/main.cpp)
target_link_libraries(FEM PRIVATE fem_core)
target_compile_options(FEM PRIVATE -Wall -Wextra)

# mesh files are opened by relative path and the csv output goes to the cwd,
# so the binary is placed in main/ and must be run from there
set_target_properties(FEM PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/main)

if(FEM_BUILD_TESTS)
	enable_testing()
	add_subdirectory(tests)
endif()
