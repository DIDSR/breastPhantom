cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(breastPhantom
    LANGUAGES CXX
)

# VTK dependencies
find_package(VTK 9.0
    COMPONENTS
        CommonComputationalGeometry
        CommonCore
        CommonDataModel
        CommonMisc
        FiltersCore
        FiltersGeneral
        ImagingHybrid
        IOXML
    QUIET
)
# for older VTK versions, pull all dependencies
if (NOT VTK_FOUND)
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
endif()
if (VTK_VERSION VERSION_LESS 6.0)
    message(FATAL_ERROR "Minimum supported version for VTK is 6")
endif()

# OpenMP
find_package(OpenMP REQUIRED)

# Optimizations and standards
set_property(GLOBAL PROPERTY CXX_STANDARD 20)
set_property(GLOBAL PROPERTY CXX_STANDARD_REQUIRED ON)
add_compile_options(-std=gnu++20)

set_property(GLOBAL PROPERTY INTERFACE_POSITION_INDEPENDENT_CODE OFF)
set_property(GLOBAL PROPERTY INTERPROCEDURAL_OPTIMIZATION true)
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)
set(CMAKE_ENABLE_EXPORTS OFF)

# Source files
include_directories("$(PROJECT_SOURCE_DIR)/src")
add_executable(${PROJECT_NAME}
    src/breastPhantom.cxx
    src/artery.cxx
    src/createArtery.cxx
    src/duct.cxx
    src/createDuct.cxx
    src/vein.cxx
    src/createVein.cxx
    src/perlinNoise.cxx
)

# Linking
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        lapack blas OpenMP::OpenMP_CXX
        ${VTK_LIBRARIES}
        boost_program_options boost_iostreams
)
if (NOT VTK_MAJOR_VERSION LESS 9)
    vtk_module_autoinit(TARGETS ${PROJECT_NAME} MODULES ${VTK_LIBRARIES})
endif()
