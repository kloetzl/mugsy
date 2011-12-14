# SEQAN Development Build System
# Platform Adapter Makefile for Linux/G++
#
# 2003/2004 by Andreas Doering
# doering@inf.fu-berlin.de
#________________________________________________________________________
#

extEXE :=

System = $(shell uname)

Compiler ?= g++
#Compiler ?= /import/gcc/bin/g++
#Compiler ?= g++-2.95
#Compiler ?= g++-3.0
#Compiler ?= g++-3.1
#Compiler ?= g++-3.2
#Compiler ?= g++-3.3
#Compiler ?= g++-3.4
#Compiler ?= g++-4.1
#Compiler ?= /import/testing/bin/g++
#Compiler ?= /bioinfo/work/data/knut/g++-4.1.2/bin/g++
#Compiler ?= /home/takifugu2/doering/usr/local/bin/g++-4.3.1

ifdef MUGSYCODEDIR
 CCFlags += -I $(MUGSYCODEDIR) -I /usr/local/projects/angiuoli/boost/include/boost-1_38 -pedantic -ftemplate-depth-200 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
else
 CCFlags += -I /usr/local/projects/angiuoli/mugsy_trunk -I /usr/local/projects/angiuoli/boost/include/boost-1_38 -pedantic -ftemplate-depth-200 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
 MUGSYCODEDIR=/usr/local/projects/angiuoli/mugsy_trunk
endif

CCFlags += -I $(MUGSYCODEDIR) -I /usr/local/projects/angiuoli/boost/include/boost-1_38 -pedantic -ftemplate-depth-200 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
#Modified by SVA
#LDFlags += /usr/local/projects/angiuoli/developer/sangiuoli/mummer/trunk/MUMmer3.20/src/tigr/*.o
LDFlags += -lstdc++

# Link against runtime library on Linux and Solaris systems
ifeq ($(System),Linux)
        LDFlags += -lrt
endif
ifeq ($(System),SunOS)
        LDFlags += -lrt
endif


# Set compiler and linker flags depending on Mode
ifeq ($(Mode),Debug)
	CCFlags += -ggdb
# Profiling -pg -W -Wall -O0
	LDFlags += -ggdb 
# Profiling	LDFlags += -pg 
endif
ifeq ($(Mode),Release)
#MODIFIED BY SVA
	CCFlags += -W -Wall -O3 
#-march=nocona -mfpmath=sse -msse2
endif
ifeq ($(Mode),Simple)
endif


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
	@mkdir -p $(1) 2> /dev/null || true
endef


# Create forwards $(1), set optional argument $(2) to "all" for rebuilding all forwards
define cmdCreateForwards
    @python misc/build_forwards.py "projects/library/seqan/$(1)" $(2)|| true
endef


# Compile $(1) to output file $(2) including files $(3)
define cmdCompile
	@echo $(Indent)compile $(1)
	@echo "  " $(Compiler) $(CCFlags) -c "$(1)" -o "$(2)" $(addprefix -I ,$(3)) 
	@$(Compiler) $(CCFlags) -c "$(1)" -o "$(2)" $(addprefix -I ,$(3))
endef


# Link $(1) to output file $(2)
define cmdLink
	@echo $(Indent)link $(2)
	@echo "  " $(Compiler) $(LDFlags) -L projects/library/apps/mugsy/ -L $(MUGSYCODEDIR)/maflib -lmaf -o "$(2)"
	#MODIFIED BY SVA
        @$(Compiler) $(1) $(LDFlags) -L projects/library/apps/mugsy/ -L $(MUGSYCODEDIR)/maflib -lmaf -o "$(2)"
endef


# Delete directories $(1) and files $(2)
define cmdClean 
	@rm -f $(2)
	@rm -r -f $(1)
endef


# Runs the application $(1)
ifeq ($(Mode),Debug)
define cmdRun 
	@valgrind -q $(1)
endef
else
define cmdRun 
	@$(1)
endef
endif


# Create documentation
define cmdDoc
	cd docs;\
	python main.py ../projects/library;\
	cd ..
endef


# Create snapshot
define cmdSnapshot
	FILENAME=Seqan_snap_`date +%Y%m%d`.zip;\
	rm -f $$FILENAME;\
	mkdir -p misc/snapshot;\
	cd projects/library/demos;\
	make clean;\
	cd ../../../misc/vs_demo_projectfile;\
	python vs_demo_projectfile.py;\
	cd ../vs_apps_projectfile;\
    python vs_apps_projectfile.py;\
	cd ../snapshot;\
	ln -s ../../projects/library/demos;\
	ln -s ../../projects/library/seqan;\
	ln -s ../../projects/library/apps;\
	ln -s ../../projects/library/cmake;\
	ln -s ../../docs/html doc;\
	echo Create $$FILENAME ...;\
	zip -qr ../../$$FILENAME * -x \*.svn\*;\
	rm demos seqan apps cmake doc;\
	cd ../..
endef

# Decompress snapshot (omit doc folder) and commit to SourceForge
define cmdSourceForge
	if [ ! -d ${constSourceForgeSVNRep} ]; then\
		echo SourceForge Repository ${constSourceForgeSVNRep} doesn\'t exist;\
		echo Please execute \"svn co https://seqan.svn.sourceforge.net/svnroot/seqan ${constSourceForgeSVNRep}\" first.;\
	else\
		FILENAME=Seqan_snap_`date +%Y%m%d`.zip;\
		if [ ! -r $$FILENAME ]; then $(call cmdSnapshot); fi;\
		unzip -qo $$FILENAME -d ${constSourceForgeSVNRep} -x doc/\*;\
		cd ${constSourceForgeSVNRep};\
		svn -q add *;\
		svn -q commit;\
	fi
endef
