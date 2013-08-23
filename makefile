CC = g++
PROCESOR = x86-64
OBJ = grid-main.o
CCOPTS = -pedantic -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fshort-enums -fno-common 
#-march=$(PROCESOR) -mtune=$(PROCESOR)
EXE = main.x
REMOVE = @rm -vf
STRIP = strip --verbose --strip-all


main:	grid-main.cpp
	$(CC) $(CCOPTS) grid-main.cpp -o main.x

$(EXE): $(OBJ) makefile 
	$(CC) $(CCOPTS) $(OBJ) -o $@ -lm
	$(STRIP) $@

*.o: *.c
	$(CC) $(CCOPTS) -c $<

clean:
	$(REMOVE) *~ *.o

veryclean: clean
	$(REMOVE) $(EXE) *.eps *.pdf *.dvi *.aux *.log *.det *.dmt *.dat
