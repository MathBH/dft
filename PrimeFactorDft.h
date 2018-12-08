#pragma once
#include <complex>
#include "PFFTMapper.h"
#include "PrimeFactors.h"
#include "DFTAlgorithm.h"
#include "BruteForceDft.h"
#include "SmallDFTs.h"

/*
	Prime Factor Fourier Transform

	Notes:
	
	An advantage of this algorithm is it allows for significantly less constraints on the sample size's properties
	This algorithm works best for sizes that split into relatively small primes and gets worse as the primes get larger

	This algorithm can be combined with other fft algorithms for computing the smaller parts who's sizes are relatively prime. This is the approach used here.
*/

class PrimeFactorDft: public DFTAlgorithm
{
private:
	PFFTMapper mapper;
	PrimeFactors primeFactorCalculator;
	DFT2 dft2;
	DFT3 dft3;
	DFT4 dft4;
	DFT5 dft5;
	DFT7 dft7;
	BruteForceDft dftBruteForce;
	std::vector<DFTAlgorithm*> subAlgorithms;

	DFTAlgorithm* getSubAlg(int size);
public:
	PrimeFactorDft();
	~PrimeFactorDft();
	std::vector<std::complex<double>> getDft(std::vector<std::complex<double>> samples) override; //TODO: refactor with object for complex list structure
};

