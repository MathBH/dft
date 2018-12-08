#pragma once
#include <vector>
#include <complex>

class DFTAlgorithm {
public:
	virtual std::vector<std::complex<double>> getDft(std::vector<std::complex<double>> samples) = 0; //TODO: refactor with object for complex list structure
};