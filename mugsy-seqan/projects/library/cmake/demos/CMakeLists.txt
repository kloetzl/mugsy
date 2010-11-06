cmake_minimum_required (VERSION 2.6) 
project (Seqan_Demos CXX)

set (CMAKE_BUILD_TYPE Debug)
#set (CMAKE_BUILD_TYPE Release)

add_subdirectory (../seqan seqan)

################################################################################
# Set path variables
################################################################################

include_directories (${SEQAN_BASE})
link_directories (${SEQAN_BASE}/lib)
if (MSVC)
	add_definitions (/W2 /wd4996 /D "_CRT_SECURE_NO_WARNINGS")
endif (MSVC)

################################################################################
# Demos
################################################################################

file(GLOB SEQAN_DEMOS ${SEQAN_BASE}/demos/[A-z]*.cpp)
foreach (DEMOFILE ${SEQAN_DEMOS})
	get_filename_component (DEMO ${DEMOFILE} NAME_WE)
	add_executable (${DEMO} ${DEMOFILE})
	if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
		target_link_libraries (${DEMO} rt)
	endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
endforeach (DEMOFILE)
