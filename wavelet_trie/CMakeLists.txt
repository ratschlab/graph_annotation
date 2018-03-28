cmake_minimum_required (VERSION 3.6.1)
project(Metannot)

enable_testing()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "DEBUG")
endif()

string(APPEND CMAKE_CXX_FLAGS " \
  -std=c++14 -Wall \
  -march=native -DNPRINT -DNDEBUG \
  -fopenmp -D_THREAD_SAFE -pthread") # -DDBGDEBUG

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
  external-libraries/sdsl-lite/include
)

#add_subdirectory()

add_executable(wtr_compress "wtr_compress.cpp" "wavelet_trie.cpp" "utils.cpp")
add_executable(wtr_merge "wtr_merge.cpp" "wavelet_trie.cpp" "utils.cpp")
add_executable(wtr_test_decompression "wtr_test_decompression.cpp" "wavelet_trie.cpp" "utils.cpp")
add_executable(pack_sd "pack_sd.cpp")
add_executable(count_unique "count_unique.cpp")
add_executable(unpack_commas "unpack_commas.cpp")

set(METALIBS
  -lsdsl -lpthread -lgmp -lboost_serialization
)

target_link_libraries(wtr_compress ${METALIBS})
target_link_libraries(wtr_merge ${METALIBS})
target_link_libraries(wtr_test_decompression ${METALIBS})
target_link_libraries(pack_sd ${METALIBS})
target_link_libraries(count_unique ${METALIBS})
target_link_libraries(unpack_commas ${METALIBS})
