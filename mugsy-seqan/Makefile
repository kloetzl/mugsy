# SEQAN Development Build System
# Makefile (Version 2)
#
# 2003/2004/2006 by Andreas Doering
# doering@inf.fu-berlin.de
#________________________________________________________________________
#
# Targets:
#	default:    			= build
# 	build: 					same as compile
#	compile:				compile and link project
#   justcompile:            same as compile but does not create forwards
#	run:					build and run
#   justrun:                run
#	clean:					delete generated files and folders
#	dumpcompilerversion:	print compiler version
#
# Arguments:
#	Platform:				name of the platform used (required!) 
#							(folders in ./platforms/)
#
#	Project:				name of project to make 
#							(folders in ./projects/*/ )
#
#	Mode:					'Debug' for debug version (debugging)
#                           'Release' for release version (fast running)
#                           'Simple' for simple version (fast compiling)
#
#	Include (opt):			a list of additional header folders or projects
#
#   Compiler (opt):         the compiler that should be used 
#                           (only for Platform = gcc)
#   
#   BuildFolder (opt):      the folder in which the project will be compiled    
#________________________________________________________________________
# Test for presence of Argument 'Platform':

# defaults
Platform ?= gcc
Mode ?= Release

#listPlatforms := $(notdir $(wildcard platforms/*))
#
#default $(MAKECMDGOALS):
#	@echo ERROR - no platform selected
#	@echo specify one of the following platforms in variable "Platform":
#	@echo $(filter-out CVS, $(listPlatforms))
#else

#________________________________________________________________________
# Define some constants:

# File extensions:
extCPP := .cpp
extOBJ := .o
extH := .h
extEXE := .exe

# Name of the library folder
constLibraryFolder := projects/library
constSourceForgeSVNRep := ../seqan.sf

#________________________________________________________________________
# Define some variables:

# Name of the build folder
BuildFolder := $(Platform)

# List of all project folders. Ignore CVS, .svn and library folders
AllProjectFolders := $(filter-out $(constLibraryFolder),$(patsubst %/.svn,,$(patsubst %/CVS,,$(wildcard projects/*))))

# List of projects
#fuer built-all hier aendern: Project ?= $(notdir $(patsubst %/CVS,,$(foreach i,$(AllProjectFolders),$(wildcard $(i)/*))))
pathProject2 := $(sort $(foreach i,$(Project),$(wildcard projects/*/$(i))))
pathProject3 := $(filter-out $(foreach i,$(Project),$(wildcard projects/library*/$(i))),$(pathProject2))
pathProject := $(pathProject3)$(foreach i,$(Project),$(wildcard projects/library/apps/$(i)))

# List of include folders
#pathInclude := $(constLibraryFolder) platforms/$(Platform)
pathInclude := $(constLibraryFolder)

#________________________________________________________________________
# Include platform dependent makefile:

platformMakefile := $(wildcard platforms/$(Platform)/platform.mak)
ifdef platformMakefile
include $(platformMakefile)
endif

#________________________________________________________________________
# Targets:

.PHONY: default build compile clean run dumpcompilerversion snapshot $(pathProject)

default: build

build: compile

dumpcompilerversion:
	$(cmdVersion)

#________________________________________________________________________
# Target compile: Compile projects

targetMakedir := $(addsuffix _makedir, $(pathProject))
targetCreateForwards := $(addsuffix _createforwards, $(Project))
targetCompileCpp := $(addsuffix _compilecpp, $(foreach i,$(pathProject),$(wildcard $(i)/*$(extCPP))))

compile: dumpcompilerversion $(pathProject)

$(pathProject): $(targetMakedir) $(targetCreateForwards) $(targetCompileCpp)
	$(call cmdLink, $(wildcard $@/$(BuildFolder)/*$(extOBJ)),$@/$(BuildFolder)/$(notdir $@)$(extEXE))

#________________________________________________________________________
# Target justcompile

targetLinkit := $(addsuffix _linkit, $(pathProject))

justcompile: dumpcompilerversion $(targetMakedir) $(targetCompileCpp) $(targetLinkit)

$(targetLinkit): %_linkit:
	$(call cmdLink, $(wildcard $*/$(BuildFolder)/*$(extOBJ)),$*/$(BuildFolder)/$(notdir $*)$(extEXE))

#________________________________________________________________________
# Create directory

$(targetMakedir): %_makedir:
	$(call cmdMkdir,$*/$(BuildFolder))

#________________________________________________________________________
# Create forwards

$(targetCreateForwards): %_createforwards:
	$(call cmdCreateForwards,$*)

#________________________________________________________________________
# Compile CPP

$(targetCompileCpp): %_compilecpp:
	$(call cmdCompile,$*,$(dir $*)$(BuildFolder)/$(patsubst %$(extCPP),%$(extOBJ),$(notdir $*)),$(pathInclude))
 
#________________________________________________________________________
# Generate HTML documentation

doc docs manual html:
	@$(call cmdDoc)

#________________________________________________________________________
# Manually Create Forwards

forwards:
	@$(call cmdCreateForwards,,all)

#________________________________________________________________________
# Create Snapshot


snapshot dist: doc
	@$(call cmdCreateForwards,,all)
	@$(call cmdSnapshot)

snapshot_nodoc:
	@$(call cmdCreateForwards,,all)
	@$(call cmdSnapshot)

#________________________________________________________________________
# Commit to SourceForge

sourceforge:
	@$(call cmdSourceForge)

#________________________________________________________________________
# Target clean: Delete generated project files
cleanPaths := $(foreach i, $(wildcard ./projects/*/*/$(BuildFolder)), "$(i)")
cleanFiles := $(foreach i, $(wildcard ./$(constLibraryFolder)/seqan/*/*_generated_forwards.h), "$(i)")

clean:
	$(call cmdClean, $(cleanPaths), $(cleanFiles))

#________________________________________________________________________
# Target run: build and run project

run: build
	@echo
	@echo =================== $(Output) START ===================
	$(call cmdRun,$(pathProject)/$(BuildFolder)/$(notdir $(pathProject))$(extEXE))
	@echo =================== $(Output)  STOP ===================
	@echo

#________________________________________________________________________

justrun:
	@echo
	@echo =================== $(Output) START ===================
	$(call cmdRun,$(pathProject)/$(BuildFolder)/$(notdir $(pathProject))$(extEXE))
	@echo =================== $(Output)  STOP ===================
	@echo


#________________________________________________________________________

#endif #of ifndef Platform
