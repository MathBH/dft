// dft.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include "WaveGenerator.h"
#include "BruteForceDft.h"
#include "PrimeFactorDft.h"
#include "PrimeFactors.h"
#include "DFTAlgorithm.h"
int main()
{
	WaveGenerator waveGen = WaveGenerator();
	BruteForceDft dftGenB = BruteForceDft();
	PrimeFactorDft pfdft = PrimeFactorDft();
	DFTAlgorithm* dftGen = &dftGenB;
	std::vector<std::complex<double>> sinwave = waveGen.sinwave(50.,100,0.005,7.);
	std::vector<std::complex<double>> coswave = waveGen.coswave(30.,100,0.01,1.);
	std::vector<std::complex<double>> sindft = dftGen->getDft(sinwave);
	std::vector<std::complex<double>> sinrebuild = dftGenB.reverseDft(sindft);
	std::vector<std::complex<double>> cosdft = dftGen->getDft(coswave);
	std::vector<std::complex<double>> cosrebuild = dftGenB.reverseDft(cosdft);
	PrimeFactors pf = PrimeFactors();
	pf.primeFactors(30);
    return 0;
}