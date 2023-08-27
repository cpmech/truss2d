cmake_minimum_required(VERSION 3.16)

project(Truss2d)

include(CTest)

set(CMAKE_CXX_STANDARD 17)

### OPTIONS ##################################################################

option(A_MAKE_VERBOSE       "Show additional messages"              OFF)
option(A_MAKE_ALL_WARNINGS  "Make with all warnings (-Wall)"        ON)
option(A_MAKE_DEBUG_SYMBOLS "Make with debug symbols (-g)"          OFF)
option(A_MAKE_PROFILING     "Make with profiling information (-pg)" OFF)
option(A_MAKE_OPTIMIZED     "Make optimized (-O3)"                  OFF)

if(A_MAKE_VERBOSE)
	set(CMAKE_VERBOSE_MAKEFILE TRUE)
endif(A_MAKE_VERBOSE)

if(A_MAKE_ALL_WARNINGS)
	add_definitions(-Wall)
endif(A_MAKE_ALL_WARNINGS)

if(A_MAKE_DEBUG_SYMBOLS)
	add_definitions(-g)
endif(A_MAKE_DEBUG_SYMBOLS)

if(A_MAKE_PROFILING)
	add_definitions(-pg)
    set(LFLAGS "${LFLAGS} -pg")
endif(A_MAKE_PROFILING)

if(A_MAKE_OPTIMIZED)
	add_definitions(-O3)
endif(A_MAKE_OPTIMIZED)

### LIBRARIES ################################################################

if(WIN32)
else(WIN32)
    # m means -lm (lib math)
	set(LIBS ${LIBS} m)
endif(WIN32)

if(WIN32)
	find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
	set(LIBS ${LIBS} MKL::MKL)
else(WIN32)
	include(FindLAPACK)
	set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
	set(LFLAGS "${LFLAGS} ${LAPACK_LINKER_FLAGS}")
endif(WIN32)

### RESULTING LIBRARY ########################################################

set(LIB_NAME "libtruss2d")

add_library(${LIB_NAME} SHARED
    src/lib/truss2d.cpp
)

target_link_libraries(${LIB_NAME} PUBLIC ${LIBS})

### SUBDIRECTORIES ###########################################################

add_subdirectory(src)
add_subdirectory(examples)
