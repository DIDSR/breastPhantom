cmake_minimum_required(VERSION 2.8)

PROJECT(breastCropProject)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories("/opt/local/include")
include_directories("/usr/local/include")
include_directories("$(PROJECT_SOURCE_DIR)")

SET(CMAKE_BUILD_TYPE "Debug")
#SET(DISABLE_IRRITATING_WARNING "-Wno-unnamed-type-template-args")
SET(CMAKE_CXX_FLAGS  "-std=c++0x ${CMAKE_CXX_FLAGS}")

LINK_DIRECTORIES("/opt/local/lib")
LINK_DIRECTORIES("/usr/local/lib")

add_executable(breastCropMain MACOSX_BUNDLE breastCrop.cxx)

target_link_libraries(breastCropMain lapack blas boost_program_options-mt ${VTK_LIBRARIES})
