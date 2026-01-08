CC=g++
CFLAGS=-g -Wall -O3 -std=c++17

SRCDIR=src
INCDIR=include
OUTDIR=obj

LOBJS = -lpthread

# Make target to build for final submission
all: pa3

pre:
	mkdir -p $(OUTDIR)

pa3: pre $(SRCDIR)/main.cpp ${OUTDIR}/matrix.o ${OUTDIR}/suraj_parser.o ${OUTDIR}/placer.o
	rm -f PA3
	$(CC) $(CFLAGS) -I$(INCDIR) $(LOBJS) $(SRCDIR)/main.cpp $(OUTDIR)/* -o PA3

${OUTDIR}/matrix.o: $(SRCDIR)/matrix.cpp $(INCDIR)/matrix.hpp
	$(CC) $(CFLAGS) -I$(INCDIR) -c $(SRCDIR)/matrix.cpp -o $(OUTDIR)/matrix.o

${OUTDIR}/suraj_parser.o: $(SRCDIR)/suraj_parser.cpp $(INCDIR)/suraj_parser.h
	$(CC) $(CFLAGS) -I$(INCDIR) -c $(SRCDIR)/suraj_parser.cpp -o $(OUTDIR)/suraj_parser.o

${OUTDIR}/placer.o: $(SRCDIR)/placer.cpp $(INCDIR)/placer.hpp
	$(CC) $(CFLAGS) -I$(INCDIR) -c $(SRCDIR)/placer.cpp -o $(OUTDIR)/placer.o

clean:
	rm -f PA3
	rm -rf $(OUTDIR)
