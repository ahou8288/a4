CC = clang
CFLAGS = -O0 -std=gnu11 -march=native -Wall -Werror
LDFLAGS = -lm -lpthread

.PHONY: all clean

all: pagerank

pagerank: pagerank.c
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

clean:
	-rm -f *.o
	-rm -f pagerank
