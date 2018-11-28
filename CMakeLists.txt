cmake_minimum_required(VERSION 2.8)

PROJECT(breastPhantomProject)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(OpenMP)
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories("$(PROJECT_SOURCE_DIR)")

add_library(perlinNoise perlinNoise.cxx)
add_library(duct duct.cxx)
add_library(artery artery.cxx)
add_library(vein vein.cxx)
add_library(createDuct createDuct.cxx)
add_library(createArtery createArtery.cxx)
add_library(createVein createVein.cxx)

SET(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_CXX_FLAGS  "-std=c++0x ${CMAKE_CXX_FLAGS}")

add_executable(breastPhantomMain breastPhantom.cxx)

target_link_libraries(breastPhantomMain perlinNoise createDuct createArtery createVein duct artery vein z lapack blas boost_program_options ${VTK_LIBRARIES})

