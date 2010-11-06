SEQAN = ..

# Link against runtime library on Linux systems
OS_NAME=$(shell uname)
ifeq ($(OS_NAME),Linux)
  LDFLAGS += -lrt
endif

CPPFLAGS += -I$(SEQAN)
CPPFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
CPPFLAGS += -O3
CPPFLAGS += -pedantic

all: $(basename $(wildcard *.cpp))
index: $(basename $(wildcard index_*.cpp))
graph: $(basename $(wildcard graph_*.cpp))
modifier: $(basename $(wildcard modifier_*.cpp))
file: $(basename $(wildcard file_*.cpp))
find: $(basename $(wildcard find_*.cpp))
clean:
	rm -f $(basename $(wildcard *.cpp))
