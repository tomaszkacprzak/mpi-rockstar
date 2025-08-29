CC=mpicxx
CFLAGS=-I../src -fopenmp
LDFLAGS=../libmpi_rockstar.a -lm -ltirpc -fopenmp

example: example.o ../libmpi_rockstar.a
	$(CC) -o $@ example.o $(LDFLAGS)

example.o: example.c
	$(CC) $(CFLAGS) -c -o $@ example.c

clean:
	rm -f example.o example