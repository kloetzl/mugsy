#
#  Possible targets:
#
#	- all (default):    builds and installs the library
#	- install:          dito
#	- test:             builds the library and the test executable
#	- clean:            removes all temporary object files and the directories
#	- cleanall:         like clean, but also removes all object files
#

#
# The seqan library location.
#

SEQAN := ../..

#
# Tools and command line arguments.
#

#MACHINE := -m32
OBJCOPY := objcopy
LDINC := $(LD) -r
GET := ./get.sh
CFLAGS := $(MACHINE) -O3 -fomit-frame-pointer -D__forceinline=__inline__ -I$(SEQAN)
CXXFLAGS := $(CFLAGS)
SYSTEM := $(shell uname)

# On Darwin use gobjcopy of MacPorts binutils
ifeq ($(SYSTEM),Darwin)
  OBJCOPY := gobjcopy
endif

ifeq ($(MACHINE),-m32)
  LDINC += -melf_i386
endif

#
# The source of the libraries.
#

VERSION_divsufsort := 1.2.3
VERSION_msufsort   := 3.1.1
URL_divsufsort     := http://homepage3.nifty.com/wpage/software/libdivsufsort-$(VERSION_divsufsort).tar.bz2
URL_msufsort       := http://www.michael-maniscalco.com/downloads/MSufSort.$(VERSION_msufsort).zip
URL_deepshallow    := http://www.mfn.unipmn.it/~manzini/lightweight/ds.tgz
URL_qsufsort       := http://www.larsson.dogma.net/suftest.c

#
# Lists of files to generate.
#

Libraries := divsufsort msufsort ds qsufsort
#Libraries := ds
Directories := $(addsuffix .d,$(Libraries))
Objects := $(addsuffix .o,$(Libraries))
WrapperCodes := $(addsuffix -api.cpp,$(Libraries))
CleanDirs := $(addsuffix -clean,$(Directories))

.PHONY: \
	clean \
	cleanall

.PRECIOUS: \
	%.tgz \
	%.tar.bz2 \
	%.zip \
	qsufsort.c \
	%.d


all: install

install: libshawarma.a
	cp $< $(SEQAN)/lib/
	
# To replace default implicit rule (%: %.o)
%: 
	
test: test.cpp libshawarma.a
	$(CXX) $(CXXFLAGS) $< libshawarma.a -o $@

# Bundle the final object file into a library.
libshawarma.a: $(Objects)
	ar rcs $@ $^

# Link all the wrappers together into one object file.
libshawarma.o: $(Objects)
	$(LDINC) $(Objects) -z muldefs -o $@



################################################################################
# divsufsort
################################################################################
divsufsort.o: divsufsort-api.o divsufsort.d/lib/divsufsort.lo
	divsufsort.d/libtool --tag=CC --mode=link $(LDINC) $^ -o $@
	$(OBJCOPY) -w -G \*seqan\* -G _ZSt\* -G _ZNSt\* $@
	
divsufsort.d/lib/divsufsort.lo: divsufsort.d
	$(MAKE) -C divsufsort.d

divsufsort-api.o: divsufsort.d

################################################################################
# MSufSort
################################################################################
msufsort.o: msufsort-api.o msufsort.d/MSufSort/MSufSort.o
	$(LDINC) $^ -o $@
	$(OBJCOPY) -w -G \*seqan\* -G _ZSt\* -G _ZNSt\* $@
	
msufsort.d/MSufSort/MSufSort.o: msufsort.d
msufsort-api.o: msufsort.d

################################################################################
# Deep-Shallow
################################################################################
ds.o: ds-api.o ds.d/ds_ssort.a
	$(LDINC) $^ -o $@
	$(OBJCOPY) -w -G \*seqan\* -G _ZSt\* -G _ZNSt\* $@
	
ds.d/ds_ssort.a: ds.d
	$(MAKE) -C ds.d ds_ssort.a CFLAGS="$(CFLAGS) -DDEBUG=0 -DNDEBUG=1"

ds-api.o: ds.d

################################################################################
# qsufsort
################################################################################
qsufsort.o: qsufsort-api.o qsufsort.d/qsufsort.o
	$(LDINC) $^ -o $@
	$(OBJCOPY) -w -G \*seqan\* -G _ZSt\* -G _ZNSt\* $@
	
qsufsort.d/qsufsort.o: qsufsort.d
qsufsort-api.o: qsufsort.d


################################################################################
# Unzip
################################################################################
divsufsort.d: divsufsort.tar.bz2
	tar xjf $<
	@mv libdivsufsort-$(VERSION_divsufsort) $@
	cd divsufsort.d; ./configure CFLAGS=$(MACHINE)

msufsort.d: msufsort.zip
	unzip $< -d $@

ds.d: ds.tgz
	mkdir $@
	tar xzfC $< $@
	
qsufsort.d: qsufsort.c
	mkdir $@
	cp $< $@/$<


################################################################################
# Download
################################################################################
divsufsort.tar.bz2:
	@$(GET) $(URL_divsufsort) $@

msufsort.zip:
	@$(GET) $(URL_msufsort) $@
	
ds.tgz:
	@$(GET) $(URL_deepshallow) $@
	
qsufsort.c:
	@$(GET) $(URL_qsufsort) $@


################################################################################
# Clean
################################################################################
clean: $(CleanDirs)
	@echo "Removing object files ..."
	@rm -f *-api.o $(Objects) libshawarma.o libshawarma.a

%-clean:
	@echo "Removing directory $(patsubst %-clean,%,$@) ..."
	@rm -rf $*

cleanall: clean
	@echo "Removing all files ..."
	@rm -f *.o *.a
