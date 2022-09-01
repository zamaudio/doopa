# This is doopa's Makefile

CXX=g++
STATIC ?= #-static -L/path/to/static/libs -lc

all:
	cd htslib && autoreconf -fi && chmod +x configure && ./configure --disable-libcurl && $(MAKE) && cd ..
	$(CXX) -g -Wall -Wno-sign-compare -O3 -funroll-loops -fomit-frame-pointer -finline-functions -std=c++11 -Ihtslib -c -o doopa.o doopa.cc
	$(CXX) $(STATIC) -o doopa doopa.o htslib/libhts.a -lz -lm -lbz2 -llzma -lpthread -lssl -lcrypto

clean:
	$(MAKE) -C htslib clean
	rm -f *.o doopa
