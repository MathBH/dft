#pragma once
#include "DFTAlgorithm.h"
#include "PrimeFactors.h"
/*
Prime Factor Fourier Transform

Implementation of PFFT algorithm from:
"Fast Fourier Transform and Convolution Algorithms" p.125-129

Notes:

An advantage of this algorithm is it allows for significantly less constraints on the sample size's properties
This algorithm works best for sizes that split into relatively small primes and gets worse as the primes get larger

This algorithm can be combined with other fft algorithms for computing the smaller parts who's sizes are relatively prime. This is the approach used here.
*/

class PFFT : public DFTAlgorithm
{
private:
	PrimeFactors pfCalulator;
public:
	PFFT();
	~PFFT();
	std::vector<std::complex<double>> getDft(std::vector<std::complex<double>> samples) override; //TODO: refactor with object for complex list structure
};

