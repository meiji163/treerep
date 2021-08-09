cmake_minimum_required(VERSION 3.14) 
project(trep)
set(CMAKE_CXX_STANDARD 11)

if(EXISTS "${CMAKE_SOURCE_DIR}/pybind11") 
	message(STATUS "Found pybind11")
	add_subdirectory(pybind11)
else()
	include(pybind11.cmake)
endif()
	
include_directories("${CMAKE_SOURCE_DIR}/src" "${CMAKE_SOURCE_DIR}/python")
#set(Python_FIND_VIRTUALENV FIRST)
#find_package(Python COMPONENTS Interpreter Development)
#set(PYTHON_EXECUTABLE "${Python_EXECUTABLE}")
#set(PYTHON_LIBRARY "${Python_LIBRARIES}")
#set(PYTHON_LIBRARY_DIR "${Python_SITELIB}")

file (GLOB SOURCE_FILES "src/*.cc")
file (GLOB HEADER_FILES "src/*.h")
file (GLOB PYTHON_FILES "python/*.cc" "python/*.h")
pybind11_add_module(trep
	${SOURCE_FILES}
	${HEADER_FILES}
	${PYTHON_FILES}
) 

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES} ${HEADER_FILES} ${PYTHON_FILES} )
#target_link_libraries(trep PUBLIC)
#install(TARGETS trep
#		COMPONENT python
#		LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}"
#	)

