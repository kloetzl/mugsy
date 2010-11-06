SEQAN = ../..

# Link against runtime library on Linux systems
OS_NAME=$(shell uname)
ifeq ($(OS_NAME),Linux)
  LDFLAGS += -lrt
endif

CPPFLAGS += -I$(SEQAN)
CPPFLAGS += -O3
CPPFLAGS += -pedantic
CPPFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

pizzachili: pizzachili.cpp $(SEQAN)/lib/libpizzachili.a
	$(CXX) $(CPPFLAGS) pizzachili.cpp -o $@ -L$(SEQAN)/lib -lpizzachili

$(SEQAN)/lib/libpizzachili.a:
	$(MAKE) -C $(SEQAN)/extra/pizzachili install
	
clean:
	rm -f $(basename $(wildcard *.cpp))
