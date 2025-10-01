# CMakeLists.txt for src directory
# This builds the extracted NPN algorithm implementation

# The extracted implementation is header-only, so we use INTERFACE library
add_library(mockturtle_extracted INTERFACE)

target_include_directories(mockturtle_extracted INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# Link with the main mockturtle library
target_link_libraries(mockturtle_extracted INTERFACE mockturtle)

# Add the extracted NPN implementation include directory
target_include_directories(mockturtle_extracted INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/algorithms/node_resynthesis
)
