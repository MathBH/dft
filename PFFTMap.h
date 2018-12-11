#pragma once
#include <vector>
#include "PrimeFactors.h"

class PFFTMap
{
public:
	virtual int getIndex(std::vector<int> indicies) = 0;
};

