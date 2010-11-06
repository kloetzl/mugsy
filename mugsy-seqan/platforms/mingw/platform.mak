# SEQAN Development Build System
# Platform Adapter Makefile for Linux/G++
#
# 2003/2004 by Andreas Doering
# doering@inf.fu-berlin.de
#________________________________________________________________________
#

extEXE :=

Compiler := mingw32-g++

CCFlags := -ftemplate-depth-200 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
LDFlags := -lstdc++


# Set compiler and linker flags depending on Mode
ifeq ($(Mode),Debug)
	CCFlags += -g -pg -W -Wall -O0
	LDFlags += -g -pg
endif
ifeq ($(Mode),Release)
	CCFlags += -W -Wall -O3
endif
ifeq ($(Mode),Simple)
endif


# Print compiler version
define cmdVersion
	@$(Compiler) -dumpversion
endef


# Create timestamp file $(1)
define cmdTouch
	@echo 1 > $(1)
endef


# Make directory $(1)
define cmdMkdir
	-@mkdir $(subst /,\,$(1))
endef


# Create forwards $(1), set optional argument $(2) to "all" for rebuilding all forwards
define cmdCreateForwards
    @misc\build_forwards.py "projects/library/seqan/$(1)" $(2)|| true
endef


# Compile $(1) to output file $(2) including files $(3)
define cmdCompile
	@echo $(Indent)compile $(1)
	@echo $(Compiler) $(CCFlags) -c "$(1)" -o "$(2)" $(addprefix -I ,$(3))
	@$(Compiler) $(CCFlags) -c "$(1)" -o "$(2)" $(addprefix -I ,$(3)) 
endef


# Link $(1) to output file $(2)
define cmdLink
	@echo $(Indent)link $(2)
	@echo $(Compiler) $(LDFlags) "$(1)" -o "$(2)"
	@$(Compiler) $(LDFlags) $(1) -o "$(2)"
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
