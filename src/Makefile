CC=gcc
BIN=$(HOME)/bin
CFLAGS=-Wall

all: $(BIN)/xyz2pov $(BIN)/pov2anim

$(BIN)/pov2anim: pov2anim.o
	$(CC) -o $(BIN)/pov2anim pov2anim.o -lm

$(BIN)/xyz2pov: xyz2pov.o atom_parser.o pov_writer.o frameCount.o
	$(CC) -o $(BIN)/xyz2pov xyz2pov.o atom_parser.o pov_writer.o frameCount.o -lm

atom_parser.o: atom_parser.c atom_parser.h

pov_writer.o: pov_writer.c pov_writer.h atom_parser.h

pov2anim.o: pov2anim.c

frameCount.o: frameCount.h frameCount.c

xyz2pov.o: pov_writer.h atom_parser.h xyz2pov.c

clean:
	rm *.o *~
