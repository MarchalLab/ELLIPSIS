cmake_minimum_required(VERSION 3.18)
project(ELLIPSIS)

set(CMAKE_CXX_STANDARD 14)

file(COPY scripts/getGeneDepth.py
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#target_compile_definitions(spliceGraph PRIVATE "BUILDDIR=\"${CMAKE_CURRENT_BINARY_DIR}\"")
add_compile_definitions("BUILDDIR=\"${CMAKE_CURRENT_BINARY_DIR}\"")

# get git sha (short) (see https://stackoverflow.com/a/3153300)
exec_program("git" "${CMAKE_CURRENT_SOURCE_DIR}"
        ARGS describe --match=NeVeRmAtCh --always --abbrev=7 --dirty
        OUTPUT_VARIABLE VERSION_SHA1 )

add_definitions("-DGIT_SHA1=\"${VERSION_SHA1}\"")

if (CMAKE_BUILD_TYPE STREQUAL Debug)
    add_definitions(-DDEBUG)
    set(CMAKE_CXX_FLAGS "-O0")
else (CMAKE_BUILD_TYPE STREQUAL Debug)
    add_definitions(-DNDEBUG)
    set(CMAKE_CXX_FLAGS "-O3")
endif (CMAKE_BUILD_TYPE STREQUAL Debug)

include_directories(${BOOST_INCLUDE_DIR})

add_subdirectory(src)