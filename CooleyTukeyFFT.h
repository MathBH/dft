#pragma once
#include "DFTAlgorithm.h"

/*
	Cooley-Tukey FFT

	Implementation of algorithm from "Numerical Recipes in C++" p.509 - p.513
	Code taken from p.513, modified to use complex<double>, std::vector<std::complex<double>>,
	and return the dft instead of replacing the data.

	TODO: book contains an efficient structure for complex number representation,
		  should consider as replacement for complex<double>
*/

class CooleyTukeyFFT : public DFTAlgorithm
{
public:
	CooleyTukeyFFT();
	~CooleyTukeyFFT();
	std::vector<std::complex<double>> getDft(std::vector<std::complex<double>> samples) override; //TODO: refactor with object for complex list structure
};