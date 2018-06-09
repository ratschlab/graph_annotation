cmake_minimum_required (VERSION 3.6.1)
cmake_policy(SET CMP0015 NEW)
set(CMAKE_VERBOSE_MAKEFILE 1)

add_library(wtr_libs STATIC
  thread_pool.cpp
  sdsl_utils.cpp
  cpp_utils.cpp
  wavelet_trie.cpp
)

link_directories(
  ../external-libraries/sdsl-lite/lib
)

include_directories(
)

target_include_directories(wtr_libs PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ../external-libraries/sdsl-lite/include
)

target_link_libraries(
  wtr_libs
  sdsl
)
