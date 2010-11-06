cmake_minimum_required (VERSION 2.6) 
project (Seqan_Tests CXX)

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
# Tests
################################################################################

get_filename_component (SEQAN_TESTS_ABS ${SEQAN_BASE}/../tests ABSOLUTE)
#if (EXISTS SEQAN_TESTS_ABS)
	file(GLOB SEQAN_TESTS ${SEQAN_TESTS_ABS}/[A-z]*)
	foreach (TEST_DIR ${SEQAN_TESTS})
		file (RELATIVE_PATH TEST_DIR_REL ${SEQAN_TESTS_ABS} ${TEST_DIR})
		file (GLOB TEST_SOURCES ${TEST_DIR}/[A-z]*.h ${TEST_DIR}/[A-z]*.cpp)
		list (SORT TEST_SOURCES)
		if (TEST_SOURCES)
			add_executable (test_${TEST_DIR_REL} ${TEST_SOURCES})
			if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
				target_link_libraries (${test_${TEST_DIR_REL}} rt)
			endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
		endif (TEST_SOURCES)
	endforeach (TEST_DIR)
#endif (EXISTS SEQAN_TESTS_ABS)
