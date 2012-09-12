CFLAGS=-std=c++0x -g $(shell pkg-config --cflags liblcb-experimental)
LDFLAGS=$(shell pkg-config --libs liblcb-experimental)

all: test_Matrix

#test: test.o
#	g++ $< -o $@ $(LDFLAGS)

#test.o: test.cpp
#	g++ -c $(CFLAGS) $< -o $@

#test_Vector: test_Vector.o
#	g++ $< -o $@ $(LDFLAGS) -lboost_system-mt -lboost_unit_test_framework-mt

#test_Vector.o: test/test_Vector.cpp
#	g++ -c $(CFLAGS) $< -o $@

test_Matrix: test_Matrix.o
	g++ $< -o $@ $(LDFLAGS) -lboost_system-mt -lboost_unit_test_framework-mt

test_Matrix.o: test_Matrix.cpp
	g++ -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o *~ a.out 

