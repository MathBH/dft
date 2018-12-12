#pragma once
#include <vector>
#include <complex>

/*
	Template for an object that calculates the dft for a given data set

	Tesla Belzile-Ha
*/

class DFTAlgorithm {
public:
	virtual std::vector<std::complex<double>> getDft(std::vector<std::complex<double>> samples) = 0; //TODO: refactor with object for complex list structure
};