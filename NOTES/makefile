mc = mpicc
flags = -g -lm -std=c99

matblocks: matblocks.c matrix.c
	$(mc) $(flags) matblocks.c matrix.c -o matblocks

clean: 
	rm -f all

