#include "stdafx.h"
#include "PFFTMapper.h"
#include <iostream>

PFFTMapper::PFFTMapper()
{
}


PFFTMapper::~PFFTMapper()
{
}

double totien(PrimeFactor pf)
{
	double value = pf.getValue();
	if (pf.getPower() == 1)
		return value - 1.;

	double base = pf.getBase();
	return value * (1. - (1. / base));
}

std::vector<int> PFFTMapper::inputMapping(int numSamples, std::vector<PrimeFactor> primeFactors)
{
	int D = primeFactors.size();
	std::vector<int> indicies = std::vector<int>();
	indicies.resize(numSamples, 0);

	int denominator = 1;
	int prime;
	int scaler;

	for (int d = 1; d < D; d++)
	{
		prime = primeFactors[d].getValue();
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

std::vector<int> PFFTMapper::outputMapping(int numSamples, std::vector<PrimeFactor> primeFactors)
{
	int D = primeFactors.size();
	std::vector<int> indicies = std::vector<int>();
	indicies.resize(numSamples, 0);

	int denominator = 1;
	PrimeFactor primeFactor;
	int primeValue;
	int scaler;
	long tot;
	long fbase;
	int ibase;
	int ti;

	int counter;

	for (int d = 1; d < D; d++)
	{
		primeFactor = primeFactors[d];
		primeValue = primeFactor.getValue();
		scaler = numSamples / primeValue;
		tot = totien(primeFactor);
		fbase = fmod(pow(scaler, tot - 1.),primeValue);
		ibase = (int)fbase;
		for (int i = 0; i < numSamples; i++)
		{
			counter = (i / denominator) % primeValue;
			ti = (ibase + (counter*primeValue))%primeValue;
			indicies[i] += counter*scaler*ti;
			indicies[i] %= numSamples;
		}
		denominator *= primeValue;
	}
	return indicies;
}
