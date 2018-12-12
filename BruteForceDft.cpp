#include "stdafx.h"
#include "BruteForceDft.h"
#define _USE_MATH_DEFINES
#include <math.h>

/*
	Note: In "Numerical Recipes in C++", the DFT algorithm does not use -i.
	Many other sources do and in testing it appeared to be more of a matter of preference
	as either way worked, only the data needs to be read differently. I have opted for negating it
	for now.

	TODO: ^Clear this up further
*/

BruteForceDft::BruteForceDft()
{

}


BruteForceDft::~BruteForceDft()
{
}


std::vector<std::complex<double>> BruteForceDft::getDft(std::vector<std::complex<double>> samples)
{
	std::complex<double> i = -1;
	i = sqrt(i);
	int numSamples = samples.size();
	std::vector<std::complex<double>> dft = std::vector<std::complex<double>>(numSamples);

	double nsth = 1. / (double)numSamples;
	for (int n = 0; n < numSamples; n++)
	{
		for (int k = 0; k < numSamples; k++)
		{
			// Note: negated i
			dft[n] += samples[k] * exp(2 * M_PI*-i*(double)k*(double)n * nsth);
		}
	}

	return dft;
}

std::vector<std::complex<double>> BruteForceDft::reverseDft(std::vector<std::complex<double>> dft)
{
	std::complex<double> i = -1;
	i = sqrt(i);
	int numSamples = dft.size();
	double numSamplesDouble = (double)numSamples;
	double nsdth = 1. / numSamplesDouble;
	std::vector<std::complex<double>> signal = std::vector<std::complex<double>>(numSamples);

	for (int n = 0; n < numSamples; n++)
	{
		for (int k = 0; k < numSamples; k++)
		{
			// Note: negated i
			signal[n] += nsdth *dft[k] * exp(-2 * M_PI*-i*(double)k*(double)n * nsdth);
		}
	}

	return signal;
}
