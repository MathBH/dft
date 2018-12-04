#pragma once
#include <vector>

class PFFTMapper
{
public:
	PFFTMapper();
	~PFFTMapper();

	std::vector<int> basicMapping(int numSamples, std::vector<int> primeBases);
};

