CC      = gcc
CFLAGS  = -g -Wall -Wno-unused-variable -O2
LIBS    = -L $(HTSDIR) -lz -lhts

# Adjust HTSDIR to point at your htslib installation #
HTSDIR = ../htslib
HTSLIB  = $(HTSDIR)/libhts.a

all: find_denovo

find_denovo: $(HTSLIB)
	$(CC) $(CFLAGS) -I $(HTSDIR) -pthread -o $@ find_denovo.c $(LIBS)


clean:
	-rm -f find_denovo

.PHONY: all clean
