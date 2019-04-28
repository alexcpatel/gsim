IDIR=include
CC=gcc
CFLAGS=-g -O3 -Wall -I$(IDIR) -std=c99
OMP=-fopenmp

ODIR=build
SDIR=src

LIBS=-lm

_DEPS = gsim.h quadtree.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = gsim.o quadtree.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BINS = gsim-omp

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS) | build_dir
	$(CC) $(OMP) -c -o $@ $< $(CFLAGS) $(LIBS)

default: all
all: clean gsim-omp

gsim-omp: $(OBJ)
	$(CC) $(OMP) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: build_dir
build_dir:
	mkdir -p $(ODIR)

.PHONY: clean
clean:
	rm -f -r $(ODIR) *~ $(BINS) $(IDIR)/*~ *.txt