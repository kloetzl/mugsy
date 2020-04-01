
.DEFAULT: all
all: libmaf.a

CPPFLAGS+=-Wall -Wextra
CFLAGS+=-O2

SOURCES=charvec.c  maf.c  multi_util.c  mz_scores.c  nib.c  seq.c  util.c
HEADERS=charvec.h  maf.h  multi_util.h  mz_scores.h  nib.h  seq.h  util.h
OBJECTS=$(SOURCES:.c=.o)

libmaf.a: $(OBJECTS)
	ar rsc $@ $^

clean:
	$(RM) -f libmaf.a $(OBJECTS)
