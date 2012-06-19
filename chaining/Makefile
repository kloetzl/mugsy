#-- Imported variables from top level makefile
# BIN_DIR AUX_BIN_DIR CXX CC CPPFLAGS LDFLAGS

CPPFLAGS = -I /usr/local/projects/angiuoli/boost/include/boost-1_38 -pedantic -ftemplate-depth-200 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O3 
#-Wall -mfpmath=sse -msse2
#-march=nocona -mfpmath=sse -msse2 
#-ggdb
# -march=pentium4
#CPPFLAGS = -ggdb -pg
#CPPFLAGS = -I /usr/local/projects/angiuoli/boost/include/boost-1_38 -pg

ifndef BIN_DIR
BIN_DIR := $(CURDIR)
endif
ifndef AUX_BIN_DIR
AUX_BIN_DIR := $(CURDIR)
endif
OBJ_RULE = $(CXX) $(CPPFLAGS) $< -static -I /usr/local/projects/angiuoli/boost/include/boost-1_38 -c -o $@
BIN_RULE = $(CXX) $(CPPFLAGS) $^ -static -L /usr/local/projects/angiuoli/boost/lib/ -lbgl-viz -lboost_graph-gcc41-mt -o $(BIN_DIR)/$@; \
           chmod 755 $(BIN_DIR)/$@
AUX_BIN_RULE = $(CXX) $(CPPFLAGS) $^ -o $(AUX_BIN_DIR)/$@; \
           chmod 755 $(AUX_BIN_DIR)/$@
VPATH := $(AUX_BIN_DIR):$(BIN_DIR)

ALL := synchain-mugsy


#-- PHONY rules --#
.PHONY: all clean


all: $(ALL)


clean:
	rm -f *.o *~
	cd $(BIN_DIR); rm -f $(ALL)
	cd $(AUX_BIN_DIR); rm -f $(ALL)




#-- not so PHONY rules --#
synchain-mugsy.o: synchain-mugsy.cpp graph.h filters.h lcbchecks.h file.h mincut.h
	$(OBJ_RULE)

synchain-mugsy: synchain-mugsy.o
	$(BIN_RULE)



