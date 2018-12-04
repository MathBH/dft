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

int main()
{
	WaveGenerator waveGen = WaveGenerator();
	BruteForceDft dftGen = BruteForceDft();
	PrimeFactorDft pfdft = PrimeFactorDft();
	std::vector<double> sinwave = waveGen.sinwave(20.,100,0.01,7.);
	std::vector<double> coswave = waveGen.coswave(30.,100,0.01,1.);
	std::vector<std::complex<double>> sindft = dftGen.getDft(sinwave);
	std::vector<std::complex<double>> sinrebuild = dftGen.reverseDft(sindft);
	std::vector<std::complex<double>> cosdft = dftGen.getDft(coswave);
	std::vector<std::complex<double>> cosrebuild = dftGen.reverseDft(cosdft);
	PrimeFactors pf = PrimeFactors();
	pf.primeFactors(30);
    return 0;
}

