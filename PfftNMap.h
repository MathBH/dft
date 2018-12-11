#pragma once
#include "PFFTMap.h"

/*
	Index mapping for n of PFFT algorithm ("Fast Fourier Transform and Convolution Algorithms", p.126)
*/

class PfftNMap : public PFFTMap
{
private:
	int N;
	int d;
	std::vector<PrimeFactor> Ni;
public:
	PfftNMap();
	PfftNMap(int size, std::vector<PrimeFactor> primeFactors);
	~PfftNMap();

	int getIndex(std::vector<int> indicies) override;
};

