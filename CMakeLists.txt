cmake_minimum_required (VERSION 3.6.1)
project(AnnoGraph)

enable_testing()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "DEBUG")
endif()

string(APPEND CMAKE_CXX_FLAGS " \
  -std=c++14 -Wall -Wextra -Werror -Wno-ignored-qualifiers \
  -msse4.2 -mavx -mavx2 -mbmi -mbmi2 -g \
  -fopenmp -D_THREAD_SAFE -pthread") # -DDBGDEBUG

# Profile build type
set(CMAKE_CXX_FLAGS_PROFILE "-pg -DNDEBUG -O2 -g")
set(CMAKE_EXE_LINKER_FLAGS_PROFILE "-pg -g")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILE "-pg -g")

# Debug build type
set(CMAKE_CXX_FLAGS_DEBUG "-O2 -g")

set(CMAKE_VERBOSE_MAKEFILE 1)
cmake_policy(SET CMP0015 NEW)

link_directories(
  external-libraries/sdsl-lite/lib
)

include_directories(
  bloom_annotator
  wavelet_trie
)

add_subdirectory(
  bloom_annotator
)

add_subdirectory(
  wavelet_trie
)

file(GLOB src_files "*.cpp")

list(FILTER src_files EXCLUDE REGEX ".*\\._.*")
list(FILTER src_files EXCLUDE REGEX ".*main.cpp")

add_library(annographlibs STATIC ${src_files})
add_executable(annograph "main.cpp")

set(METALIBS
  -lhts -lz -lbz2
  -lssl -lcrypto -llzma
  -lsdsl -lgmp
)

if(BUILD_STATIC)
  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} -static)
  set(METALIBS ${METALIBS} -lcurl)
endif()

#target_include_directories(
#  annographlibs PUBLIC
#)

target_link_libraries(
  annographlibs bloom_annotator wtr_libs ${METALIBS}
)

target_link_libraries(
  annograph annographlibs ${METALIBS}
)


#-------------------
# Unit Tests
#-------------------
# Download and unpack googletest at configure time
if(NOT EXISTS ${CMAKE_BINARY_DIR}/googletest-download)
  configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
  endif()
endif()

# Prevent overriding the parent project's compiler/linker settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

add_subdirectory(
  ${CMAKE_BINARY_DIR}/googletest-src
  ${CMAKE_BINARY_DIR}/googletest-build
  EXCLUDE_FROM_ALL
)

file(GLOB unit_tests_files "tests/*.cpp")
list(FILTER unit_tests_files EXCLUDE REGEX ".*\\._.*")

if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

include_directories(
  .
)

link_directories(
  ${CMAKE_BINARY_DIR}
)

add_executable(unit_tests ${unit_tests_files})

target_link_libraries(unit_tests gtest_main gtest annographlibs ${METALIBS})

add_test(NAME unit_tests COMMAND unit_tests)
