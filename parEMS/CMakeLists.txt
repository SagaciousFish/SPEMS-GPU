cmake_minimum_required(VERSION 3.0)
cmake_policy(VERSION 3.0)

project(ParallelEMS)

SET(VERSION_MAJOR "0")
SET(VERSION_MINOR "1")
SET(VERSION_PATCH "0")

if(CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -pthread -fopenmp -ltbb") 
endif(CMAKE_COMPILER_IS_GNUCXX)


##########################
message("##########")
option(PROFILE "PROFILE" OFF)
if(PROFILE)
	add_definitions(-DPROFILE)
	message(STATUS "Allow PROFILE")
else()
	message(STATUS "Not Allow PROFILE")
endif()

option(BREAKDOWN "BREAKDOWN" ON)
if(BREAKDOWN )
	add_definitions(-DBREAKDOWN)
	message(STATUS "Allow BREAKDOWN")
else()
	message(STATUS "Not Allow BREAKDOWN")
endif()


option(SEQINTERSECTION "SEQINTERSECTION" ON)
if(SEQINTERSECTION )
	add_definitions(-DSEQINTERSECTION)
	message(STATUS "Allow Seq intersection")
else()
	message(STATUS "Not Allow Seq intersection")
endif()
##########################


include_directories (./include)
include_directories (./include/ems2)

# It will find all files with specific extension names
file(GLOB SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

add_library(${PROJECT_NAME} ${SRC})

add_subdirectory(test)

