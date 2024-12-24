cmake_minimum_required(VERSION 3.30)
cmake_policy(VERSION 3.30)
include(FetchContent)

project(ParallelEMS VERSION 0.1.0 LANGUAGES CXX)

add_library(${PROJECT_NAME})
file(GLOB SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
target_sources(${PROJECT_NAME} PRIVATE ${SRC})

message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
FetchContent_Declare(
        TBB
        GIT_REPOSITORY https://github.com/uxlfoundation/oneTBB.git
        GIT_TAG v2022.0.0
        GIT_SHALLOW ON
)
option(TBB_TEST OFF)
FetchContent_MakeAvailable(TBB)
target_include_directories(${PROJECT_NAME} PUBLIC BEFORE ${TBB_SOURCE_DIR})
#find_package(TBB REQUIRED)

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)
find_package(OpenMP)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -O3)

# compile in C++ 20
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)
set(CMAKE_CXX_STANDARD 20)

target_link_libraries(${PROJECT_NAME} PUBLIC TBB::tbb Threads::Threads)
if (OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif ()


##########################
message("##########")
option(PROFILE "PROFILE" OFF)
if (PROFILE)
    add_definitions(-DPROFILE)
    message(STATUS "Allow PROFILE")
else ()
    message(STATUS "Not Allow PROFILE")
endif ()

option(BREAKDOWN "BREAKDOWN" ON)
if (BREAKDOWN)
    add_definitions(-DBREAKDOWN)
    message(STATUS "Allow BREAKDOWN")
else ()
    message(STATUS "Not Allow BREAKDOWN")
endif ()


option(SEQINTERSECTION "SEQINTERSECTION" ON)
if (SEQINTERSECTION)
    add_definitions(-DSEQINTERSECTION)
    message(STATUS "Allow Seq intersection")
else ()
    message(STATUS "Not Allow Seq intersection")
endif ()
##########################

include_directories(./include)
include_directories(./include/ems2)

# We adding a lot of tests here
add_subdirectory(test)

