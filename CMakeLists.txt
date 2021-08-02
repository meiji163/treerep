cmake_minimum_required(VERSION 2.8.12)
project(trep)
set(CMAKE_CXX_STANDARD 11)
add_subdirectory(pybind11)
include_directories("${CMAKE_SOURCE_DIR}/src" "${CMAKE_SOURCE_DIR}/python")

file (GLOB SOURCE_FILES "src/*.cc")
file (GLOB HEADER_FILES "src/*.h")
file (GLOB PYTHON_FILES "python/*.cc" "python/*.h")
pybind11_add_module(trep
	${SOURCE_FILES}
	${HEADER_FILES}
	${PYTHON_FILES}
) 

