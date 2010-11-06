# SEQAN Development Build System
# Platform Adapter Makefile for Linux/G++ + LiSA
#
# 2003/2004 by Andreas Doering
# doering@inf.fu-berlin.de
#________________________________________________________________________
#

LISA_DIR := /home/takifugu/mbauer/research/lisa
LISA_LIB_DIR := ${LISA_DIR}/lib
LISA_INC_DIR := ${LISA_DIR}/incl/

LEDA_DIR := /import/leda/LEDA-5.0-complete-FC2_386-g++-3.3.3-std
LEDA_LIB_DIR := ${LEDA_DIR}
LEDA_INC_DIR := ${LEDA_DIR}/incl_new

VIENNA_DIR := /home/arabidopsis/lisa/software/ViennaRNA-1.6
VIENNA_LIB_DIR := ${VIENNA_DIR}/lib
VIENNA_INC_DIR := ${VIENNA_DIR}/include/ViennaRNA

extEXE := 

Compiler := g++
#Compiler := /import/gcc/bin/g++
#Compiler := g++-2.95
#Compiler := g++-3.0
#Compiler := g++-3.2
#Compiler := g++-3.3
#Compiler := g++-3.4
#Compiler := g++-4.1
#Compiler := /import/testing/bin/g++
#Compiler := /bioinfo/work/data/knut/g++-4.1.2/bin/g++

# Print compiler version
define cmdVersion
	@$(Compiler) -dumpversion
endef

# Create timestamp file $(1)
define cmdTouch
	-@touch $(1)
endef

# Make directory $(1)
define cmdMkdir
	@mkdir $(1) 2> /dev/null || true
endef

# Create forwards $(1), set optional argument $(2) to "all" for rebuilding all forwards
define cmdCreateForwards
    @python misc/build_forwards.py "projects/library/seqan/$(1)" $(2)|| true
endef


ifeq ($(Mode),Debug)
#	CCFlags := -ggdb -O0 -Wall -pedantic -ftemplate-depth-200
CCFlags := -Wall -pedantic -O0 -g -pg -ftemplate-depth-200

	LDFlags := -g -pg
else
	CCFlags := -Wall -pedantic -O3 -ftemplate-depth-200
	LDFlags :=
endif

# Compile $(1) to output file $(2) including files $(3)
define cmdCompile
	@echo $(Indent)compile $(1)
	@echo "  " $(Compiler) $(CCFlags) "$(1)" -o "$(2)" $(addprefix -I ,$(3) ) -I ${LISA_INC_DIR} -I ${LEDA_INC_DIR} -I ${VIENNA_INC_DIR}
	@$(Compiler) $(CCFlags) -c "$(1)" -o "$(2)" $(addprefix -I ,$(3) ) -I ${LISA_INC_DIR} -I ${LEDA_INC_DIR} -I ${VIENNA_INC_DIR}
endef

# Link $(1) to output file $(2)
define cmdLink
	@echo $(Indent)link $(2)
	@echo "  " $(Compiler) $(LDFlags) $(1) -o "$(2)" -L ${LISA_LIB_DIR} -L ${LEDA_LIB_DIR} -L ${VIENNA_LIB_DIR} -lstdc++ -lrt -lLISALEDA -lLISABASE -lLISA -lcommon -lRNA -lP -lG -lL -lm

	@$(Compiler) $(LDFlags) $(1) -o "$(2)" -L ${LISA_LIB_DIR} -L ${LEDA_LIB_DIR} -L ${VIENNA_LIB_DIR} -lstdc++ -lrt -lLISALEDA -lLISABASE -lLISA -lLISAVienna -lcommon -lRNA -lP -lG -lL -lm
endef

# Delete directories $(1) and files $(2)
define cmdClean 
	@rm -f $(2)
	@rm -r -f $(1)
endef

# Runs the application $(1)
define cmdRun 
	@$(1)
endef
