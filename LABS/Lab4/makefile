mc = mpicc
flags = -g -lm -std=c99

all: l4

l4: main.o matrix.o
	$(mc) $(flags) main.o matrix.o -o l4

main.o: main.c
	$(mc) -c $(flags) main.c

matrix.o: matrix.c
	$(mc) -c $(flags) matrix.c

clean: 
	rm -f *.o
	rm l4
