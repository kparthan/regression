CFLAGS=-std=c++0x $(shell pkg-config --cflags liblcb-0.1.0)
LDFLAGS=$(shell pkg-config --libs liblcb-0.1.0)

all: regression test_Vector test_Matrix

regression: regression.o
	g++ $< -o $@ $(LDFLAGS)

regression.o: regression.cpp
	g++ -c $(CFLAGS) $< -o $@

test_Vector: test_Vector.o
	g++ $< -o $@ $(LDFLAGS) -lboost_system-mt -lboost_unit_test_framework-mt

test_Vector.o: test_Vector.cpp
	g++ -c $(CFLAGS) $< -o $@

test_Matrix: test_Matrix.o
	g++ $< -o $@ $(LDFLAGS) -lboost_system-mt -lboost_unit_test_framework-mt

test_Matrix.o: test_Matrix.cpp
	g++ -c $(CFLAGS) $< -o $@

