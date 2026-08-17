# Catch2 is fetched at configure time rather than vendored, so the repository
# holds none of its source. v3 is the current major version and is a compiled
# library rather than the single header of v2.
include(FetchContent)

FetchContent_Declare(
	Catch2
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG v3.5.2
	GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(Catch2)

add_executable(fem_tests
	test_quadrature.cpp
	test_basis.cpp
	test_convergence.cpp
)

# Catch2WithMain supplies main(), so the test files contain only test cases
target_link_libraries(fem_tests PRIVATE fem_core Catch2::Catch2WithMain)
target_compile_options(fem_tests PRIVATE -Wall -Wextra)

# the convergence tests read meshes by relative path, exactly as the solver does
set_target_properties(fem_tests PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/main)

# registers each TEST_CASE with ctest individually
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(fem_tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/main)
