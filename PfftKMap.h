#pragma once
#include "PFFTMap.h"
#include "LCSolver.h"

/*
	Index mapping for k of PFFT algorithm ("Fast Fourier Transform and Convolution Algorithms", p.126)
*/

class PfftKMap : public PFFTMap
{
private:
	int N;
	int d;
	std::vector<PrimeFactor> Ni;
	LCSolver lcSolver;
public:
	PfftKMap();
	PfftKMap(int size, std::vector<PrimeFactor> primeFactors);
	~PfftKMap();

	int getIndex(std::vector<int> indicies) override;
};

