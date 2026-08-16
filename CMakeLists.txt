cmake_minimum_required(VERSION 3.16)
project(FEMwithCpp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Eigen3 REQUIRED NO_MODULE)

file(GLOB FEM_SOURCES ${CMAKE_SOURCE_DIR}/source/*.cpp)

add_executable(FEM ${FEM_SOURCES} ${CMAKE_SOURCE_DIR}/main/main.cpp)

target_include_directories(FEM PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(FEM PRIVATE Eigen3::Eigen)

# the stricter set the .vscode config used adds ~400 more warnings, almost all
# sign-conversion noise: -Wpedantic -Wshadow -Wformat=2 -Wcast-align
# -Wconversion -Wsign-conversion -Wnull-dereference
target_compile_options(FEM PRIVATE -Wall -Wextra)

# mesh files are opened by relative path and CSVs are written to the CWD,
# so the binary is placed in main/ and must be run from there
set_target_properties(FEM PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/main)
