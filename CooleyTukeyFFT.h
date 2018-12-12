#pragma once
#include "DFTAlgorithm.h"
/*
	CooleyTukeyFFT.h
	Edited from "Numerical Recipes in C++"

	Tesla Belzile-Ha
*/

/*
	Class calculating the dft for a sample set using the Cooley-Tukey algorithm.
	
	Edits:
		- Uses complex<double> and std::vector<std::complex<double>> instead of double and arrays
		- Returns a seperately built dft instead of replacing the data.

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