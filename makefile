VER=2.0
DIST=qubic$(VER)
PROGS=qubic
SRCS=struct.c read_array.c make_graph.c get_options.c fib.c write_block.c cluster.c main.c expand.c
OBJS=$(SRCS:.c=.o) 
CC=gcc

LDFLAGS= -lm -L/usr/local/gsl/latest/lib -lgsl -lgslcblas
CFLAGS=-O3 -Wall -ansi -I/usr/local/gsl/latest/include  -DVER=$(VER)

all: $(PROGS)

${PROGS}: $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

.o:
	$(CC) $(CFLAGS) $< -o $@ 

clean:
	rm -f $(PROGS)
	rm -f *.o
	rm -f data/*.rules
	rm -f data/*.chars
	rm -f data/*.blocks
	rm -f data/*.expansion

dist:
	$(MAKE) clean
	cd .. && tar czvf $(DIST).tar.gz $(DIST)/

test: 
	$(MAKE)
	./${PROGS} -i data/example 
