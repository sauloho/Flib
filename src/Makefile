#################################################################
#       Makefile for Flib                                       #
#################################################################

#### CFLAGS     	Flags to pass to C compiler.
CFLAGS          	= -O2


all:			Flib extract filterlib filterlib2 Flib_Enrich parse_hhr get_length LibValidator

clean:
			/bin/rm -rf *.o Flib Flib_Enrich LibValidator get_length filterlib filterlib2 parse_hhr extract

install: 
			/bin/cp Flib extract filterlib filterlib2 Flib_Enrich parse_hhr get_length LibValidator ../bin

Flib:			slink.h jacobi.h Flib.o coordcol_FS2.o super4.o jacobi.o sparse.o
			$(CC) $(CFLAGS) Flib.o coordcol_FS2.o super4.o jacobi.o sparse.o -lm -o Flib

Flib_Enrich:		slink.h jacobi.h main_FS4.o coordcol_FS2.o super4.o jacobi.o
			$(CC) main_FS4.o coordcol_FS2.o super4.o jacobi.o -lm -o Flib_Enrich -Wall -O2

LibValidator:		slink.h jacobi.h coordcol_FS2.o super4.o LibValidator.o 
			$(CC) $(CFLAGS) LibValidator.o coordcol_FS2.o super4.o jacobi.o sparse.o -lm -o LibValidator

get_length:		get_length.c
			$(CC) $(CFLAGS) get_length.c -o get_length

LibValidator.o:		slink.h LibValidator.c
			$(CC) $(CFLAGS) -c LibValidator.c

filterlib:		filter_lib.c
			$(CC) $(CFLAGS) filter_lib.c -o filterlib

filterlib2:		filter_lib2.c
			$(CC) $(CFLAGS) filter_lib2.c -o filterlib2

parse_hhr:		parse_hhr.c
			$(CC) $(CFLAGS) parse_hhr.c -o parse_hhr

extract:		slink.h coordcol_FS3.o extract.o
			$(CC) $(CFLAGS) coordcol_FS3.o extract.o -lm -o extract

extract.o:		slink.h extract.c
			$(CC) $(CFLAGS) -c extract.c

coordcol_FS3.o:		slink.h coordcol_FS3.c
			$(CC) $(CFLAGS) -c coordcol_FS3.c

flib.o:			slink.h Flib.c
			$(CC) $(CFLAGS) -c Flib.c

coordcol_FS2.o:		slink.h coordcol_FS2.c
			$(CC) $(CFLAGS) -c coordcol_FS2.c

sparse.o:		sparse.c
			$(CC) $(CFLAGS) -lm -c sparse.c

super4.o:		slink.h jacobi.h super4.c
			$(CC) $(CFLAGS) -c super4.c

jacobi.o:		jacobi.h jacobi.c
			$(CC) $(CFLAGS) -c jacobi.c

main_FS4.o:		slink.h main_FS4.c
			$(CC) -c main_FS4.c -Wall -O2

