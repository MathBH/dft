#pragma once
#include <vector>
#include <complex>
#include "DFTAlgorithm.h"

/*
	Brute Force DFT
	Tesla Belzile-Ha

	A brute force implementation based on equations from based on equations
	from "Numerical Recipes in C++" p.508.

*/

class BruteForceDft: public DFTAlgorithm
{
public:
	BruteForceDft();
	~BruteForceDft();
	std::vector<std::complex<double>> getDft(std::vector<std::complex<double>> samples) override; //TODO: refactor with object for complex list structure
	std::vector<std::complex<double>> reverseDft(std::vector<std::complex<double>> dft); //TODO: refactor with object for complex list structure
};

