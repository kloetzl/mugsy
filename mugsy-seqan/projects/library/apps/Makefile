SEQAN_BASE = ..

# Link against runtime library on Linux systems
OS_NAME=$(shell uname)
ifeq ($(OS_NAME),Linux)
  LDFLAGS += -lrt
endif

tbb_root?=../extra/tbb

#check, if tbb_root is not absolute path (the filter keeps only /* paths)
ifeq ($(filter /% $(SLASH)%, $(subst :, ,$(tbb_root)) ),)
    # also changes related variables like work_dir
    override tbb_root := $(CWD)$(SLASH)..
    export TBB21_INSTALL_DIR := $(tbb_root)
endif

# explicitly compile for a 32 or 64 bit platform
#CXXFLAGS += -m32
#CXXFLAGS += -m64

CXXFLAGS += -I$(SEQAN_BASE)
CXXFLAGS += -O9 -march=nocona -mfpmath=sse -msse2
#CXXFLAGS += -O0 -g
CXXFLAGS += -pedantic -W -Wall
CXXFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 ${CXXEXTRAFLAGS}

TARGETS = dfi/dfi seqan_tcoffee/seqan_tcoffee seqcons/seqcons razers/paramChooser razers/razers pair_align/pair_align micro_razers/micro_razers tree_recon/tree_recon mugsy/mugsy
mugsy: mugsy/mugsy

all:           check_seqan_base $(TARGETS)
dfi:           check_seqan_base dfi/dfi 
razers:        check_seqan_base razers/razers razers/paramChooser
micro_razers:  check_seqan_base micro_razers/micro_razers
seqan_tcoffee: check_seqan_base seqan_tcoffee/seqan_tcoffee 
seqcons:       check_seqan_base seqcons/seqcons
pair_align:    check_seqan_base pair_align/pair_align
tree_recon:    check_seqan_base tree_recon/tree_recon

check_seqan_base:
	@if [ ! -d "$(SEQAN_BASE)/seqan" ]; then \
		echo "The directory $(SEQAN_BASE)/seqan could not be found!"; \
		exit 1; \
	fi

clean:
	rm -f $(TARGETS) $(TARGETS:=.o)
