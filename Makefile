CC = gcc
CFLAGS = -Wall -c -std=c11 -O3
OBJECTS = glla.o
RED   = '\033[0;31m'
GREEN = '\033[0;32m'
NOCOL = '\033[0m'

all: $(OBJECTS) glla.h Makefile

.depend:
	gcc -M *.c > .depend

install: all
	ar -cvq libglla.a $(OBJECTS)
	mv libglla.a /usr/local/lib/
	cp glla.h /usr/local/include/

test: glla_test.c glla.c glla.h
	$(CC) glla_test.c glla.o -o glla_test
	@if ./glla_test; then echo -e "${GREEN}Tests succeeded.${NOCOL}"; else echo -e "${RED}Tests failed.${NOCOL}"; fi

clean:
	rm $(OBJECTS) glla_test

include .depend
