add_library(bloom_annotator STATIC
  dbg_bloom_annotator.cpp
  hashers.cpp
  ../serialization.cpp
)

include_directories(
  ../external-libraries/rollinghashcpp
)

target_include_directories(bloom_annotator PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(
  bloom_annotator
)
