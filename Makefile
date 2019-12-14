# This is doopa's Makefile

CXX=g++

all:
	$(MAKE) -C htslib
	$(CXX) -g -Wall -Wno-sign-compare -O3 -funroll-loops -fomit-frame-pointer -finline-functions -fopenmp -std=c++11 -Ihtslib -c -o doopa.o doopa.cc
	$(CXX) -o doopa doopa.o htslib/libhts.a -fopenmp -lz -lm -lbz2 -llzma -lpthread -lcurl

clean:
	rm *.o doopa
