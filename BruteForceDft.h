#pragma once
#include <vector>
#include <complex>

class BruteForceDft
{
public:
	BruteForceDft();
	~BruteForceDft();
	std::vector<std::complex<double>> getDft(std::vector<double> samples); //TODO: refactor with object for complex list structure
	std::vector<std::complex<double>> reverseDft(std::vector<std::complex<double>> dft); //TODO: refactor with object for complex list structure
};

