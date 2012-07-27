CFLAGS=-std=c++0x -g $(shell pkg-config --cflags liblcb-experimental)
LDFLAGS=$(shell pkg-config --libs liblcb-experimental)

all: regression test_Vector test_Matrix

regression: regression.o
	g++ $< -o $@ $(LDFLAGS)

regression.o: regression.cpp Message.h
	g++ -c $(CFLAGS) $< -o $@

test_Vector: test_Vector.o
	g++ $< -o $@ $(LDFLAGS) -lboost_system-mt -lboost_unit_test_framework-mt

test_Vector.o: test_Vector.cpp
	g++ -c $(CFLAGS) $< -o $@

test_Matrix: test_Matrix.o
	g++ $< -o $@ $(LDFLAGS) -lboost_system-mt -lboost_unit_test_framework-mt

test_Matrix.o: test_Matrix.cpp
	g++ -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o regression a.out test_Vector test_Matrix

