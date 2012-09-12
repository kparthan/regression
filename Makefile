CFLAGS=-std=c++0x -g $(shell pkg-config --cflags liblcb-experimental)
LDFLAGS=$(shell pkg-config --libs liblcb-experimental)

all: regression 

regression: regression.o
	g++ $< -o $@ $(LDFLAGS)

regression.o: regression.cpp *.h
	g++ -c $(CFLAGS) $< -o $@

#test_Vector: test_Vector.o
#	g++ $< -o $@ $(LDFLAGS) -lboost_system-mt -lboost_unit_test_framework-mt

#test_Vector.o: test/test_Vector.cpp
#	g++ -c $(CFLAGS) $< -o $@

#test_Matrix: test_Matrix.o
#	g++ $< -o $@ $(LDFLAGS) -lboost_system-mt -lboost_unit_test_framework-mt

#test_Matrix.o: test/test_Matrix.cpp
#	g++ -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o *~ regression a.out test/test_Vector test/test_Matrix *.png temp/* *.p

