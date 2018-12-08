#pragma once
#include <vector>
#include <complex>
#include "DFTAlgorithm.h"

class BruteForceDft: public DFTAlgorithm
{
public:
	BruteForceDft();
	~BruteForceDft();
	std::vector<std::complex<double>> getDft(std::vector<std::complex<double>> samples) override; //TODO: refactor with object for complex list structure
	std::vector<std::complex<double>> reverseDft(std::vector<std::complex<double>> dft); //TODO: refactor with object for complex list structure
};

