CC      = g++
CFLAGS  = -g
RM      = rm -f


default: all

all: build run

build: main.cpp OTL.cpp OTL.h
	$(CC) $(CFLAGS) -o main.o main.cpp OTL.cpp

run: main.o
	@./main.o testargument

clean veryclean:
	$(RM) main.o
