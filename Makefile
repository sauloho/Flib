all: Flib extract filterlib filterlib2 Flib_Enrich parse_hhr get_length LibValidator

CFLAGS=-O2

Flib: slink.h jacobi.h Flib.o coordcol_FS2.o super4.o jacobi.o sparse.o
	cc $(CFLAGS) Flib.o coordcol_FS2.o super4.o jacobi.o sparse.o -lm -o Flib

Flib_Enrich: slink.h jacobi.h main_FS4.o coordcol_FS2.o super4.o jacobi.o
	cc main_FS4.o coordcol_FS2.o super4.o jacobi.o -lm -o Flib_Enrich -Wall -O2

LibValidator: slink.h jacobi.h coordcol_FS2.o super4.o LibValidator.o 
	cc $(CFLAGS) LibValidator.o coordcol_FS2.o super4.o jacobi.o sparse.o -lm -o LibValidator

get_length: get_length.c
	cc $(CFLAGS) get_length.c -o get_length

LibValidator.o: slink.h LibValidator.c
	cc $(CFLAGS) -c LibValidator.c

filterlib: filter_lib.c
	cc $(CFLAGS) filter_lib.c -o filterlib

filterlib2: filter_lib2.c
	cc $(CFLAGS) filter_lib2.c -o filterlib2

parse_hhr: parse_hhr.c
	cc $(CFLAGS) parse_hhr.c -o parse_hhr

extract: slink.h coordcol_FS3.o extract.o
	cc $(CFLAGS) coordcol_FS3.o extract.o -lm -o extract
extract.o: slink.h extract.c
	cc $(CFLAGS) -c extract.c
coordcol_FS3.o:slink.h coordcol_FS3.c
	cc $(CFLAGS) -c coordcol_FS3.c
flib.o: slink.h Flib.c
	cc $(CFLAGS) -c Flib.c
coordcol_FS2.o:slink.h coordcol_FS2.c
	cc $(CFLAGS) -c coordcol_FS2.c
sparse.o: sparse.c
	cc $(CFLAGS) -lm -c sparse.c
super4.o:slink.h jacobi.h super4.c
	cc $(CFLAGS) -c super4.c
jacobi.o: jacobi.h jacobi.c
	cc $(CFLAGS) -c jacobi.c

main_FS4.o: slink.h main_FS4.c
	cc -c main_FS4.c -Wall -O2

clean:
	rm -rf *.o Flib parse_hhr filterlib filterlib2 Flib_Enrich extract get_length
