IDIR=include
CC=gcc
CFLAGS=-Wall -I$(IDIR)
OMP=-fopenmp

ODIR=build
SDIR=src

LIBS=-lm

_DEPS = gsim.h monitor.h cycletimer.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = gsim.o monitor.o cycletimer.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BINS = gsim-bad

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS) | build_dir
	$(CC) $(OMP) -c -o $@ $< $(CFLAGS) $(LIBS)

default: all
all: clean gsim-bad

gsim-bad: $(OBJ)
	$(CC) $(OMP) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: build_dir
build_dir:
	mkdir -p $(ODIR)

.PHONY: clean
clean:
	rm -f -r $(ODIR) *~ $(BINS) $(IDIR)/*~ *.txt