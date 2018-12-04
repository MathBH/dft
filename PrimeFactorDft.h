#pragma once
#include <complex>
#include "PFFTMapper.h"
#include "PrimeFactors.h"

/*
	Prime Factor Fourier Transform

	Notes:
	
	An advantage of this algorithm is it allows for significantly less constraints on the sample size's properties
	This algorithm works best for sizes that split into relatively small primes and gets worse as the primes get larger

	This algorithm can be combined with other fft algorithms for computing the smaller parts who's sizes are relatively prime. This is the approach used here.
*/

class PrimeFactorDft
{
private:
	PFFTMapper mapper;
	PrimeFactors primeFactorCalculator;
public:
	PrimeFactorDft();
	~PrimeFactorDft();
	std::vector<std::complex<double>> getDft(std::vector<double> samples); //TODO: refactor with object for complex list structure
};

