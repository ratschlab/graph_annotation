cmake_minimum_required (VERSION 3.6.1)
project(Metannot)

enable_testing()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "DEBUG")
endif()

string(APPEND CMAKE_CXX_FLAGS " \
  -std=c++14 -Wall \
  -march=native -DNPRINT \
  -fopenmp -D_THREAD_SAFE -pthread -g") # -DDBGDEBUG

# Profile build type
set(CMAKE_CXX_FLAGS_PROFILE "-pg -DNDEBUG -O2 -g")
set(CMAKE_EXE_LINKER_FLAGS_PROFILE "-pg -g")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILE "-pg -g")

# Debug build type
set(CMAKE_CXX_FLAGS_DEBUG "-O3 -g")

set(CMAKE_VERBOSE_MAKEFILE 1)
cmake_policy(SET CMP0015 NEW)

link_directories(
  external-libraries/sdsl-lite/lib
)

include_directories(
)

set(SDSLINCLUDE "external-libraries/sdsl-lite/include")

#add_subdirectory()

add_library(wtr_libs STATIC "thread_pool.cpp")
#add_executable(wtr_compress "wtr_compress.cpp" "wavelet_trie.cpp")
#add_executable(wtr_merge "wtr_merge.cpp" "wavelet_trie.cpp")
#add_executable(wtr_test_decompression "wtr_test_decompression.cpp" "wavelet_trie.cpp")
add_executable(pack_sd "pack_sd.cpp")
add_executable(count_unique "count_unique.cpp")
#add_executable(unpack_commas "unpack_commas.cpp")

set(METALIBS
  -lsdsl -lpthread -lgmp -lboost_serialization
)

#add_library(wavelet_trie STATIC external-libraries/sdsl-lite/lib/libsdsl.a "wtr_compress.cpp" "wavelet_trie.cpp" "thread_pool.cpp" "wtr_merge.cpp" "wtr_test_decompression.cpp")
#target_include_directories(wtr_compress PUBLIC ${SDSLINCLUDE})
#target_include_directories(wtr_merge PUBLIC ${SDSLINCLUDE})
#target_include_directories(wtr_test_decompression PUBLIC ${SDSLINCLUDE})
target_include_directories(pack_sd PUBLIC ${SDSLINCLUDE})
target_include_directories(count_unique PUBLIC ${SDSLINCLUDE})
#target_include_directories(unpack_commas PUBLIC ${SDSLINCLUDE})
#target_link_libraries(wtr_compress wtr_libs ${METALIBS})
#target_link_libraries(wtr_merge wtr_libs ${METALIBS})
#target_link_libraries(wtr_test_decompression wtr_libs ${METALIBS})
target_link_libraries(pack_sd ${METALIBS})
target_link_libraries(count_unique ${METALIBS})
#target_link_libraries(unpack_commas ${METALIBS})

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
file(GLOB wavelet_trie_files "wavelet_trie.cpp")
list(FILTER unit_tests_files EXCLUDE REGEX ".*\\._.*")

if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

include_directories(
  .
  ${SDSLINCLUDE}
)

link_directories(
  ${CMAKE_BINARY_DIR}
)

add_executable(unit_tests ${unit_tests_files} ${wavelet_trie_files})

target_link_libraries(unit_tests gtest_main gtest wtr_libs ${METALIBS})

add_test(NAME unit_tests COMMAND unit_tests)
