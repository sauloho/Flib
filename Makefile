# Compiler Options
CFLAGS=-O2

# Folders
BIN=bin
SRC=src

all: Flib extract filterlib filterlib2 Flib_Enrich parse_hhr get_length LibValidator

Flib: $(SRC)/slink.h $(SRC)/jacobi.h $(SRC)/Flib.o $(SRC)/coordcol_FS2.o $(SRC)/super4.o $(SRC)/jacobi.o $(SRC)/sparse.o
	$(CC) $(CFLAGS) $(SRC)/Flib.o $(SRC)/coordcol_FS2.o $(SRC)/super4.o $(SRC)/jacobi.o $(SRC)/sparse.o -lm -o $(BIN)/Flib

Flib_Enrich: $(SRC)/slink.h $(SRC)/jacobi.h $(SRC)/main_FS4.o $(SRC)/coordcol_FS2.o $(SRC)/super4.o $(SRC)/jacobi.o
	$(CC) $(SRC)/main_FS4.o $(SRC)/coordcol_FS2.o $(SRC)/super4.o $(SRC)/jacobi.o -lm -o $(BIN)/Flib_Enrich -Wall -O2

LibValidator: $(SRC)/slink.h $(SRC)/jacobi.h $(SRC)/coordcol_FS2.o $(SRC)/super4.o $(SRC)/LibValidator.o 
	$(CC) $(CFLAGS) $(SRC)/LibValidator.o $(SRC)/coordcol_FS2.o $(SRC)/super4.o $(SRC)/jacobi.o $(SRC)/sparse.o -lm -o $(BIN)/LibValidator

get_length: $(SRC)/get_length.c
	$(CC) $(CFLAGS) $(SRC)/get_length.c -o $(BIN)/get_length

LibValidator.o: $(SRC)/slink.h $(SRC)/LibValidator.c
	$(CC) $(CFLAGS) -c $(SRC)/LibValidator.c

filterlib: $(SRC)/filter_lib.c
	$(CC) $(CFLAGS) $(SRC)/filter_lib.c -o $(BIN)/filterlib

filterlib2: $(SRC)/filter_lib2.c
	$(CC) $(CFLAGS) $(SRC)/filter_lib2.c -o $(BIN)/filterlib2

parse_hhr: $(SRC)/parse_hhr.c
	$(CC) $(CFLAGS) $(SRC)/parse_hhr.c -o $(BIN)/parse_hhr

extract: $(SRC)/slink.h $(SRC)/coordcol_FS3.o $(SRC)/extract.o
	$(CC) $(CFLAGS) $(SRC)/coordcol_FS3.o $(SRC)/extract.o -lm -o $(BIN)/extract
extract.o: $(SRC)/slink.h $(SRC)/extract.c
	$(CC) $(CFLAGS) -c $(SRC)/extract.c
coordcol_FS3.o: $(SRC)/slink.h $(SRC)/coordcol_FS3.c
	$(CC) $(CFLAGS) -c $(SRC)/coordcol_FS3.c
flib.o: $(SRC)/slink.h $(SRC)/Flib.c
	$(CC) $(CFLAGS) -c $(SRC)/Flib.c
coordcol_FS2.o: $(SRC)/slink.h $(SRC)/coordcol_FS2.c
	$(CC) $(CFLAGS) -c $(SRC)/coordcol_FS2.c
sparse.o: $(SRC)/sparse.c
	$(CC) $(CFLAGS) -lm -c $(SRC)/sparse.c
super4.o: $(SRC)/slink.h $(SRC)/jacobi.h $(SRC)/super4.c
	$(CC) $(CFLAGS) -c $(SRC)/super4.c
jacobi.o: $(SRC)/jacobi.h $(SRC)/jacobi.c
	$(CC) $(CFLAGS) -c $(SRC)/jacobi.c

main_FS4.o: $(SRC)/slink.h $(SRC)/main_FS4.c
	$(CC) -c $(SRC)/main_FS4.c -Wall -O2

clean:
	rm -rf $(SRC)/*.o $(BIN)/Flib $(BIN)/Flib_Enrich $(BIN)/LibValidator $(BIN)/get_length $(BIN)/filterlib $(BIN)/filterlib2 $(BIN)/parse_hhr $(BIN)/extract
