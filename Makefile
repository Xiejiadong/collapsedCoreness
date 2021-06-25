CC=g++ -march=native -O3
CFLAGS=-c -I. -std=c++11 -Wfatal-errors

all: core

core: main.o heap.o treap.o glist.o
	$(CC) main.o heap.o treap.o glist.o -o core
	rm *.o

main.o: main.cc
	$(CC) $(CFLAGS) main.cc -o main.o

heap.o: heap.cc
	$(CC) $(CFLAGS) heap.cc -o heap.o

treap.o: treap.cc
	$(CC) $(CFLAGS) treap.cc -o treap.o

glist.o: glist.cc
	$(CC) $(CFLAGS) glist.cc -o glist.o

