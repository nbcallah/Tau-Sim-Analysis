CC=gcc
CXX=g++
OPT=-O3
CFLAGS=$(OPT)
CPPFLAGS=$(OPT) -std=c++11 -fopenmp
CPPLFLAGS=-fopenmp

all: chisq_spectrum_fit

debug: OPT=-g
debug: all

chisq_spectrum_fit: xorshift.o chisq_spectrum_fit.o
	$(CXX) $(CPPLFLAGS) -o chisq_spectrum_fit chisq_spectrum_fit.o xorshift.o

chisq_spectrum_fit.o: chisq_spectrum_fit.cpp
	$(CXX) $(CPPFLAGS) -c -o chisq_spectrum_fit.o chisq_spectrum_fit.cpp

xorshift.o: xorshift.c xorshift.h
	$(CC) $(CFLAGS) -c -o xorshift.o xorshift.c

clean:
	rm -rf xorshift.o chisq_spectrum_fit.o chisq_spectrum_fit