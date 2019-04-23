IDIR=include
CC=gcc
CFLAGS=-g -O3 -Wall -I$(IDIR) -std=c99

ODIR=build
SDIR=src

LIBS=-lm

_DEPS = gsim.h body.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = gsim.o body.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BINS = gsim

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS) | build_dir
	$(CC) -c -o $@ $< $(CFLAGS)

default: all
all: clean gsim

gsim: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: build_dir
build_dir:
	mkdir -p $(ODIR)

.PHONY: clean
clean:
	rm -f -r $(ODIR) *~ $(BINS) $(IDIR)/*~ *.txt