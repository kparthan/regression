CFLAGS=-std=c++0x -g $(shell pkg-config --cflags liblcb-experimental)
LDFLAGS=$(shell pkg-config --libs liblcb-experimental)

all: regression 

regression: regression.o
	g++ $< -o $@ $(LDFLAGS)

regression.o: regression.cpp *.h
	g++ -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o *~ regression a.out test/test_Vector test/test_Matrix *.eps *.p temp/*.eps temp/*.p

