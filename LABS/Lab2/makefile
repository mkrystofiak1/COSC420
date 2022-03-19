mc = mpicc
flags = -g -lm -std=c99

main: main.o matrix.o
	$(mc) $(flags) main.o matrix.o -o main

main.o: main.c
	$(mc) -c $(flags) main.c

matrix.o: matrix.c
	$(mc) -c $(flags) matrix.c

clean: 
	rm -f *.o
	rm main

