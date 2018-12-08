CXX=g++

all:
	$(MAKE) -C htslib
	$(CXX) -g -Wall -O2 -Ihtslib -c -o doopa.o doopa.cc
	$(CXX) -o doopa doopa.o htslib/libhts.a -lz -lm -lbz2 -llzma -lpthread

clean:
	rm *.o doopa
