# Modified for library-use by Hiroki Shibata
# external/Re-Pair/CMakeLists.txt
cmake_minimum_required(VERSION 3.16)
project(RePair LANGUAGES CXX)

add_library(re_pair STATIC
        rp_lib.cpp
)

target_compile_features(re_pair PUBLIC cxx_std_17)

target_include_directories(re_pair
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/internal>
        $<INSTALL_INTERFACE:include/RePair>
)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

target_compile_options(re_pair PRIVATE
        "$<$<CONFIG:Release>:-Ofast;-fstrict-aliasing;-DNDEBUG>"
        "$<$<CONFIG:RelWithDebInfo>:-Ofast;-fstrict-aliasing>"
        "$<$<CONFIG:Debug>:-O0;-g>"
)

add_library(RePair::re_pair ALIAS re_pair)
