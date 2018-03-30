CC=gcc
CXX=g++
OPT=-O3
OMP=-fopenmp
CFLAGS=$(OPT)
CPPFLAGS=$(OPT) -std=c++11 $(OMP) `root-config --cflags`
CPPLFLAGS=$(OMP) -lgsl -lgslcblas `root-config --libs`

all: chisq_spectrum_fit arrival_time_fixed_eff pse_sum_fit

debug: OPT=-g
debug: all

single: OMP=
single: all

pse_sum_fit: xorshift.o pse_sum_fit.o
	$(CXX) -o pse_sum_fit pse_sum_fit.o xorshift.o $(CPPLFLAGS)
    
pse_sum_fit.o: pse_sum_fit.cpp
	$(CXX) $(CPPFLAGS) -c -o pse_sum_fit.o pse_sum_fit.cpp

chisq_spectrum_fit: xorshift.o chisq_spectrum_fit.o
	$(CXX) -o chisq_spectrum_fit chisq_spectrum_fit.o xorshift.o $(CPPLFLAGS)

chisq_spectrum_fit.o: chisq_spectrum_fit.cpp
	$(CXX) $(CPPFLAGS) -c -o chisq_spectrum_fit.o chisq_spectrum_fit.cpp

xorshift.o: xorshift.c xorshift.h
	$(CC) $(CFLAGS) -c -o xorshift.o xorshift.c
    
arrival_time_fixed_eff: arrival_time_fixed_eff.cpp
	$(CXX) $(CPPFLAGS) -o arrival_time_fixed_eff arrival_time_fixed_eff.cpp $(CPPLFLAGS)

clean:
	rm -rf xorshift.o chisq_spectrum_fit.o chisq_spectrum_fit arrival_time_fixed_eff