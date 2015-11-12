cmake_minimum_required(VERSION 2.8)

PROJECT(breastMassProject)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(OpenMP)
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#include_directories("/opt/local/include")
#include_directories("/usr/local/include")
include_directories("$(PROJECT_SOURCE_DIR)")

SET(CMAKE_BUILD_TYPE "Debug")
#SET(DISABLE_IRRITATING_WARNING "-Wno-unnamed-type-template-args")
SET(CMAKE_CXX_FLAGS  "-std=c++0x ${CMAKE_CXX_FLAGS}")

#LINK_DIRECTORIES("/opt/local/lib")

#add_executable(breastMassMain MACOSX_BUNDLE breastMass.cxx)
add_executable(breastMassMain breastMass.cxx)

target_link_libraries(breastMassMain lapack blas boost_program_options ${VTK_LIBRARIES})
#target_link_libraries(breastMassMain PoissonReconstruction lapack blas boost_program_options-mt ${VTK_LIBRARIES})
