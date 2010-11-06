cmake_minimum_required (VERSION 2.6) 
project (Seqan_All CXX)

set (CMAKE_BUILD_TYPE Debug)
#set (CMAKE_BUILD_TYPE Release)

add_subdirectory (apps)
add_subdirectory (demos)
add_subdirectory (tests)
