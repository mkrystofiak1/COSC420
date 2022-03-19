mc = mpicc
flags = -g -lm -std=c99

pagerank: pagerank.o matrix.o
	$(mc) $(flags) pagerank.o matrix.o -o pagerank

testCondense: testCondense.o matrix.o
	$(mc) $(flags) testCondense.o matrix.o -o testCondense

write: write.o
	$(mc) $(flags) write.o -o write

pagerank.o: pagerank.c
	$(mc) -c $(flags) pagerank.c

testCondense.o: testCondense.c
	$(mc) -c $(flags) testCondense.c

write.o: write.c
	$(mc) -c $(flags) write.c

matrix.o: matrix.c
	$(mc) -c $(flags) matrix.c

run: pagerank.c 
	mpirun -n 3 pagerank

clean:
	rm -f *.o
	rm -f pagerank
	rm -f testCondense
	rm -f write
