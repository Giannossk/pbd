cmake_minimum_required(VERSION 3.16)
project(pbd VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Dependency: LinearMath
include(FetchContent)

option(BUILD_TESTING "Build tests" OFF)

if(NOT DEFINED FETCHCONTENT_SOURCE_DIR_LINEARMATH AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../LinearMath/CMakeLists.txt")
    set(FETCHCONTENT_SOURCE_DIR_LINEARMATH "${CMAKE_CURRENT_SOURCE_DIR}/../LinearMath" CACHE PATH "Source directory for LinearMath" FORCE)
endif()

FetchContent_Declare(
    LinearMath
    GIT_REPOSITORY https://github.com/Giannossk/LinearMath.git
    GIT_TAG        main
    GIT_SHALLOW    FALSE
)

# Ensure CMake checks for remote updates when configuring
set(FETCHCONTENT_UPDATES_DISCONNECTED_LINEARMATH OFF CACHE BOOL "Always check for the latest commit on configure" FORCE)

FetchContent_MakeAvailable(LinearMath)

# Collect source files
file(GLOB_RECURSE PBD_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
)

# PBD library
add_library(pbd ${PBD_SOURCES})
add_library(pbd::pbd ALIAS pbd)

target_include_directories(pbd
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(pbd
    PUBLIC
        LinearMath
)

if(MSVC)
    target_compile_options(pbd PRIVATE /EHsc /W4)
else()
    target_compile_options(pbd PRIVATE -Wall -Wextra)
endif()

# Testing
option(PBD_BUILD_TESTING "Build PBD tests" ON)

if(PBD_BUILD_TESTING)
    enable_testing()
    add_subdirectory(testing)
endif()
