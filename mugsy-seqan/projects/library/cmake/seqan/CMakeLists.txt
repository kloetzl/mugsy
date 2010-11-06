cmake_minimum_required (VERSION 2.6) 
project (Seqan CXX)

#set_property (GLOBAL PROPERTY ALLOW_DUPLICATE_CUSTOM_TARGETS TRUE)

################################################################################
# Set path variables
################################################################################

if (NOT SEQAN_BASE)
	set (SEQAN_BASE ../.. PARENT_SCOPE)
endif (NOT SEQAN_BASE)

set (SEQAN_TARGET "SEQAN")
#set_property (GLOBAL PROPERTY ALLOW_DUPLICATE_CUSTOM_TARGETS TRUE)

################################################################################
# Library sources
################################################################################

#if (NOT TARGET ${SEQAN_TARGET})
#get_property(IS_DEFINED GLOBAL PROPERTY SEQAN_TARGET_DEFINED)
#if (NOT ${IS_DEFINED})
#	set_property(GLOBAL PROPERTY SEQAN_TARGET_DEFINED TRUE)
	
	# Find all library headers
	get_filename_component (SEQAN_BASE_ABS ${SEQAN_BASE}/seqan ABSOLUTE)
	file (GLOB SEQAN_HEADERS ${SEQAN_BASE_ABS}/[A-z]*/[A-z]*.h)
	list (SORT SEQAN_HEADERS)
#	add_custom_target (${SEQAN_TARGET} SOURCES ${SEQAN_HEADERS})
	
	# Group library headers into modules
	foreach (SEQAN_HEADER ${SEQAN_HEADERS})
		file (RELATIVE_PATH SEQAN_HEADER_REL ${SEQAN_BASE_ABS} ${SEQAN_HEADER})
		get_filename_component (SEQAN_MODULE ${SEQAN_HEADER_REL} PATH)
		source_group (${SEQAN_MODULE} FILES ${SEQAN_HEADER})
	endforeach (SEQAN_HEADER)

#endif (NOT TARGET ${SEQAN_TARGET})
#endif (NOT ${IS_DEFINED})

################################################################################
# Set compiler flags
################################################################################

# Make sure the compiler can find include files from the SeqAn library. 
include_directories (${SEQAN_BASE})

if (CMAKE_COMPILER_IS_GNUCXX)
	set (CMAKE_CXX_WARNING_LEVEL 4 PARENT_SCOPE)
endif (CMAKE_COMPILER_IS_GNUCXX)

if (MSVC)
	add_definitions (/W2 /wd4996 /D "_CRT_SECURE_NO_WARNINGS")
endif (MSVC)
