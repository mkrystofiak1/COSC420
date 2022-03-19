all: hello prime

hello: hello.c
	mpicc -std=c99 -o hello hello.c

prime: prime.c
	mpicc -std=c99 -lm -o prime prime.c

clean:
	rm -f *.o all


