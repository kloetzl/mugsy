#
#  Possible targets:
#
#	- all (default):    builds and installs the index
#	- install:          dito
#	- test:             builds the index and the test executable
#	- clean:            removes all temporary binary files and the directories
#	- cleanall:         like clean, but also removes all binary files
#

#
# Tools and command line arguments.
#

OBJCOPY := objcopy
LDINC := $(LD) -r
GET := ./get.sh
CPPFLAGS := -Wall -pedantic
SYSTEM := $(shell uname)

# On Darwin use gobjcopy of MacPorts binutils
ifeq ($(SYSTEM),Darwin)
  OBJCOPY := gobjcopy
endif

#
# The source of the libraries.
#

Host := http://pizzachili.dcc.uchile.cl
Rel := indexes

#
# The seqan library location.
#

SEQAN := ../..

#
# Pizza & Chili libraries that go into the final bundle.
#

AF := Alphabet-Friendly_FM-Index/af-index.tgz

# Not a reusable library.
CSA := Compact_Suffix_Array/CSA.tgz

CCSA := Compressed_Compact_Suffix_Array/CCSA.tgz

FM := FM-indexV2/fmindexV2.tgz

# Buggy: redefines standard system name ulong from
# /usr/include/sys/types.h:151
LZ := LZ-index/LZ-index-4.tar.gz

RSA := Repair_Suffix_Array/RPSA.tgz

RLFM := Run_Length_FM_index/RLFM.tgz

SA := Suffix_Array/SAc.tgz

SADA := Compressed_Suffix_Array/sada_csa.tgz

# Use version 1 instead of 2 (the latter makes problems).
SSA := Succinct_Suffix_Array/SSA.tgz

#
# Define all C++ libraries here like this:
#
# MyLib_type := CPP
#
# C libraries can be omitted, or may be defined as "C" instead of "CPP".
#

AF_type := CPP
CSA_type := C
CCSA_type := CPP
FM_type := C
LZ_type := C
RSA_type := CPP
RLFM_type := CPP
SA_type := C
SADA_type := C
SSA_type := CPP

#
# Lists of files to generate.
#

# SSA and RLFM make problems.
#Libraries := AF CCSA FM RSA RLFM SA SADA SSA
Libraries := AF CCSA FM RSA SA SADA
LibLang := $(addsuffix _type,$(Libraries))
Directories := $(addsuffix .d,$(Libraries))
Tarballs := $(addsuffix .tgz,$(Libraries))
Objects := $(addsuffix .o,$(Libraries))
WrapperCodes := $(addprefix pizzachili_api_,$(addsuffix .cpp,$(Libraries)))
WrapperObjects := $(addsuffix -wrap.o,$(Libraries))
CleanDirs := $(addsuffix -clean,$(Directories))

.PHONY: \
	clean \
	cleanall

.PRECIOUS: \
	%.tgz \
	%.d

################################################################################
#
# Steps necessary to compile the libraries:
#
#  1. Download the index library.
#  2. Rename the downloaded archive to %.tgz where % is a "friendly" name.
#  3. Create a subdirectory and move the archive there.
#  4. Untar the archive.
#  5. Remove the %.tgz file.
#  6. Build the source by recursively calling make.
#  7. Extract all the content of .a to a temporary directory.
#  8. Compile the wrapper.
#  9. Re-link the wrapper.o with all the .o files from the temporary directory.
# 10. Mark all the private Pizza & Chili symbols in wrapper.o local.
# 11. Link all the wrappers together into one object file.
# 12. Bundle the final object file into a library.
#
################################################################################

all: install

install: libpizzachili.a
	cp $< $(SEQAN)/lib/
	
# To replace default implicit rule (%: %.o)
%: 
	
test: test.cpp libpizzachili.a
	$(CXX) $(CPPFLAGS) -I$(SEQAN) $< libpizzachili.a -o $@

# 12. Bundle the final object file into a library.
libpizzachili.a: libpizzachili.o
	ar rcs $@ $<

# 11. Link all the wrappers together into one object file.
libpizzachili.o: $(WrapperObjects)
	$(LDINC) $(WrapperObjects) -z muldefs -o $@

#  8. Compile the wrapper.
#  9. Re-link the wrapper.o with all the .o files from the temporary directory.
# 10. Mark all the private Pizza & Chili symbols in wrapper.o local.
%-wrap.o: %.o pizzachili_api.cpp interface.h
	$(CXX) $(CPPFLAGS) -I$(SEQAN) -c -o $*-tmp.o \
		-DPCINDEX_NAME=PizzaChiliApi$* \
		-DPCINDEX_$* -DPCINDEX_$($*_type) pizzachili_api.cpp
	$(LDINC) -x $*-tmp.o $*.d/tmp/*.o -o $@
	$(OBJCOPY) -w -G \*PizzaChili\* -G _ZSt\* -G _ZNSt\* $@
	
#. 6. Build the source by recursively calling make.
#  7. Extract all the content of .a to a temporary directory.
%.o: %.d
	@$(MAKE) -C $< clean
	@$(MAKE) -C $<
	@cd $< && \
		mkdir -p tmp && \
		mv *.a tmp && \
		cd tmp && \
		ar x *.a && \
		rm *.a && \
		cd ../.. # Not strictly necessary in this particular version of make.

FM.o: FM.d
	@$(MAKE) -C FM.d/fmindexV2 clean
	@$(MAKE) -C FM.d/fmindexV2
	@cd FM.d && \
		mkdir -p tmp && \
		mv fmindexV2/fm_index.a tmp && \
		cd tmp && \
		ar x *.a && \
		rm *.a && \
		cd ../.. # Not strictly necessary in this particular version of make.

#. 3. Create a subdirectory and move the archive there.
#. 4. Untar the archive.
%.d: %.tgz
	@mkdir -p $@
	@cd $@; tar xzf ../$*.tgz

#. 1. Download the index library.
#. 2. Rename the downloaded archive to %.tgz where % is a "friendly" name.
%.tgz:
	@$(GET) $(Host)/$(Rel)/$($*)
	@-mv $(notdir $($*)) $@

clean: $(CleanDirs)
	@echo "Removing object files ..."
	@rm -f *-tmp.o *-wrap.o libpizzachili.o libpizzachili.a

%-clean:
	@echo "Removing directory $(patsubst %-clean,%,$@) ..."
	@rm -rf $*

cleanall: clean
	@echo "Removing all files ..."
	@rm -f *.o *.a
