
# Makefile for pdbcat

# need a C++ compilier
CC = /usr/bin/g++
CFLAGS = -c -O

# if your compilier complains about strcasecmp (gcc did)
#   uncomment the next line
#DEF = -DNOSTRCASECMP

pdbcat: pdbcat.o PDBData.o Common.o
	$(CC) pdbcat.o PDBData.o Common.o -o pdbcat

pdbcat.o: pdbcat.C
	$(CC) $(CFLAGS) pdbcat.C

PDBData.o: PDBData.C PDBData.h
	$(CC) $(CFLAGS) PDBData.C

Common.o: Common.C Common.h
	$(CC) $(CFLAGS) Common.C $(DEF)

clean:
	rm -f *.o pdbcat
