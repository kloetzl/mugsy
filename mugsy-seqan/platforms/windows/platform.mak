# SEQAN Development Build System
# Platform Adapter Makefile for Visual C++
#
# 2003/2004 by Andreas Doering
# doering@inf.fu-berlin.de
#________________________________________________________________________
#

extEXE := .exe
SHELL := $(ComSpec)

# Print compiler version
define cmdVersion
endef

# Create timestamp file $(1)
define cmdTouch
	@echo 1 > $(1)
endef

# Make directory $(1)
define cmdMkdir
	-@mkdir $(subst /,\,$(1))
endef

# Create forwards $(1)
#define cmdCreateForwards
#	@misc\build_forwards.py "projects/library/seqan/$(1)"
#endef

define cmdCreateForwards
endef

# Detect VS Version
ifndef Version
	ifdef VS71COMNTOOLS
		Version=7
	endif
	ifdef VS80COMNTOOLS
		Version=8
	endif
	ifdef VS90COMNTOOLS
		Version=9
	endif
endif

ifeq ($(Version),7)
	VSCOMNTOOLS := $(VS71COMNTOOLS)
	
	ifeq ($(Mode),Debug)
		CompileFlags = /Op- /EHsc /D "DEBUG" /D "WIN32" /Zi /GR /W2 /Zc:wchar_t
		LinkFlags = /debug
	endif
	ifeq ($(Mode),Release)
		  # CompileFlags = /Op- /EHsc /Ox
		CompileFlags = /Op- /EHsc /D "NDEBUG" /D "WIN32" /O2 /Ob2 /W2 /Zc:wchar_t
		  #LinkFlags = /debug 
		  LinkFlags = 
	endif
	ifeq ($(Mode),Simple)
		CompileFlags = /Op- /EHsc /D "NDEBUG" /D "WIN32" /W2 /Zc:wchar_t
		LinkFlags =
	endif
else
	ifeq ($(Version),8)
		VSCOMNTOOLS := $(VS80COMNTOOLS)
	endif
	ifeq ($(Version),9)
		VSCOMNTOOLS := $(VS90COMNTOOLS)
	endif
	
	ifeq ($(Mode),Debug)
		CompileFlags = /EHsc /MTd /Zi /D "_CONSOLE" /D "_CRT_SECURE_NO_DEPRECATE" /D "DEBUG" /D "WIN32" /Zc:wchar_t /W2 /wd4996
		  # CompileFlags = /Za /Op- /EHsc /Zi /GR
		LinkFlags = /DEBUG
	endif
	ifeq ($(Mode),Release)
		CompileFlags = /EHsc /Ox /D "NDEBUG" /D "_CRT_SECURE_NO_DEPRECATE" /D "WIN32" /Zc:wchar_t /W2 /wd4996
		  # CompileFlags = /Za /Op- /EHsc /O2 /Zi /Ob2 /D "NDEBUG"
		LinkFlags = 
	endif
	ifeq ($(Mode),Simple)
		CompileFlags = /EHsc /D "NDEBUG" /D "_CONSOLE" /D "_CRT_SECURE_NO_DEPRECATE" /D "WIN32" /Zc:wchar_t /W2 /wd4996
		LinkFlags =
	endif
endif

# Compile $(1) to output file $(2) including files $(3)
define cmdCompile
	@echo $(Indent)compile $(1)
	@"$(VSCOMNTOOLS)\vsvars32.bat" && CL /c $(addprefix /I, $(3)) $(addprefix -I ,"projects/libraryDiff") /Fo$(2) $(1) $(CompileFlags)
endef

# Link $(1) to output file $(2)
define cmdLink
	@echo $(Indent)link $(2)
	@"$(VSCOMNTOOLS)\vsvars32.bat" && LINK "/OUT:$(2)" $(1) $(LinkFlags)
endef

# Delete directories $(1) and files $(2)
define cmdClean 
    -@del /Q $(subst /,\,$(2))
	-@rmdir /Q /S $(subst /,\,$(1))
endef

# Runs the application $(1)
define cmdRun 
	@$(subst /,\,$(1))
endef

# Create documentation
define cmdDoc
	@cd docs
	@main.py ../projects/library;\
	@cd ..
endef
