# Makefile

2mm: 2mm.c 2mm.h
	gcc -o 2mm 2mm.c

clean:
	rm -rf *.o
