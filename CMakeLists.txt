cmake_minimum_required(VERSION 3.16)
project(Project)

set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --show-leak-kinds=definite,indirect,possible")
set(MEMORYCHECK_SUPPRESSIONS_FILE "${PROJECT_SOURCE_DIR}/zscripts/valgrind.supp" CACHE FILEPATH "" FORCE)
include(CTest)

set(CMAKE_CXX_STANDARD 20)

### options ##################################################################

option(A1_OPTIMIZED "Make optimized (-O3)" OFF)

if(A1_OPTIMIZED)
    add_definitions(-O3)
else()
    add_definitions(-g -Wall)
endif()

# dependencies ###############################################################################

find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
include(zscripts/FindLACLIB.cmake)
include_directories(${LACLIB_INCS})

### library ##################################################################

add_library(truss2d SHARED src/lib/truss2d.cpp)
target_compile_definitions(truss2d PUBLIC USE_MKL)
target_link_libraries(truss2d PUBLIC MKL::MKL ${LACLIB_LIBS})

### SUBDIRECTORIES ###########################################################

add_subdirectory(src)
add_subdirectory(examples)
