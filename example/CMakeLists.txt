cmake_minimum_required(VERSION 2.8)

PROJECT(simpleExample_Project)

set(CMAKE_BUILD_TYPE Release)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -Werror -O2")
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(${OpenMP_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src/ )
file(GLOB SimpleCT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../src/*)
add_library(lSimpleCT STATIC ${SimpleCT_SRC})

SET(example_SRC main.cpp)

ADD_EXECUTABLE(simpleExample ${example_SRC})
	
TARGET_LINK_LIBRARIES(simpleExample
	rt
	mpi
	lSimpleCT
	)

INSTALL(TARGETS simpleExample DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})