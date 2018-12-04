#include "stdafx.h"
#include "PFFTMapper.h"
#include <iostream>

PFFTMapper::PFFTMapper()
{
}


PFFTMapper::~PFFTMapper()
{
}

std::vector<int> PFFTMapper::basicMapping(int numSamples, std::vector<int> primeBases)
{
	int D = primeBases.size();
	std::vector<int> indicies = std::vector<int>();
	indicies.resize(numSamples, 0);

	int denominator = 1;
	int prime;
	int scaler;

	for (int d = 0; d < D; d++)
	{
		prime = primeBases[d];
		scaler = numSamples / prime;
		for (int i = 0; i < numSamples; i++)
		{
			indicies[i] += ((i / denominator) % prime)*scaler;
			indicies[i] %= numSamples;
		}
		denominator *= prime;
	}

	return indicies;
}
