#pragma once
#include <vector>
#include "PrimeFactors.h"

class PFFTMapper
{
public:
	PFFTMapper();
	~PFFTMapper();

	std::vector<int> basicMapping(int numSamples, std::vector<PrimeFactor> primeBases);
};

