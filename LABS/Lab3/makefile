mc = mpicc
flags = -g -lm -std=c99

all: l1 l3

l1: main.o matrix.o
	$(mc) $(flags) main.o matrix.o -o l1

l3: newmain.o matrix.o
	$(mc) $(flags) newmain.o matrix.o -o l3

main.o: main.c
	$(mc) -c $(flags) main.c

newmain.o: newmain.c
	$(mc) -c $(flags) newmain.c

matrix.o: matrix.c
	$(mc) -c $(flags) matrix.c

clean: 
	rm -f *.o
	rm l1
	rm l3
