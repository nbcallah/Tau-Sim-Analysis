CC=gcc
CXX=g++
CFLAGS=-O3
CPPFLAGS=-O3 -std=c++11

all: chisq_spectrum_fit

chisq_spectrum_fit: xorshift.o chisq_spectrum_fit.o
	$(CXX) -o chisq_spectrum_fit chisq_spectrum_fit.o xorshift.o

chisq_spectrum_fit.o: chisq_spectrum_fit.cpp
	$(CXX) $(CPPFLAGS) -c -o chisq_spectrum_fit.o chisq_spectrum_fit.cpp

xorshift.o: xorshift.c xorshift.h
	$(CC) $(CFLAGS) -c -o xorshift.o xorshift.c