#pragma once
#include <vector>
#include "PrimeFactors.h"

class PFFTMapper
{
public:
	PFFTMapper();
	~PFFTMapper();

	std::vector<int> inputMapping(int numSamples, std::vector<PrimeFactor> primeFactors);
	std::vector<int> outputMapping(int numSamples, std::vector<PrimeFactor> primeFactors);
};

