cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

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
if (VTK_VERSION VERSION_LESS 5.10)
    message(FATAL_ERROR "Minimum supported version for VTK is 5.10")
endif()

# OpenMP
find_package(OpenMP REQUIRED)

# Optimizations and standards
set_property(GLOBAL PROPERTY CXX_STANDARD 14)
set_property(GLOBAL PROPERTY CXX_STANDARD_REQUIRED ON)
add_compile_options(-std=gnu++14)

set_property(GLOBAL PROPERTY INTERFACE_POSITION_INDEPENDENT_CODE OFF)
set_property(GLOBAL PROPERTY INTERPROCEDURAL_OPTIMIZATION true)
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# Source files
include_directories("$(PROJECT_SOURCE_DIR)")
add_executable(${PROJECT_NAME}
    breastPhantom.cxx
    artery.cxx
    createArtery.cxx
    duct.cxx
    createDuct.cxx
    vein.cxx
    createVein.cxx
    perlinNoise.cxx
)

# Linking
target_link_libraries(${PROJECT_NAME}
    PRIVATE z lapack blas boost_program_options ${VTK_LIBRARIES} ${OpenMP_CXX_FLAGS}
)
if (NOT VTK_MAJOR_VERSION LESS 9)
    vtk_module_autoinit(TARGETS ${PROJECT_NAME} MODULES ${VTK_LIBRARIES})
endif()
