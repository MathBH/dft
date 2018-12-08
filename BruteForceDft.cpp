#include "stdafx.h"
#include "BruteForceDft.h"
#define _USE_MATH_DEFINES
#include <math.h>

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

	for (int n = 0; n < numSamples; n++)
	{
		for (int k = 0; k < numSamples; k++)
		{
			dft[n] += samples[k] * exp(2 * M_PI*i*(double)k*(double)n / (double)numSamples);
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
	std::vector<std::complex<double>> signal = std::vector<std::complex<double>>(numSamples);

	for (int n = 0; n < numSamples; n++)
	{
		for (int k = 0; k < numSamples; k++)
		{
			signal[n] += (1./numSamplesDouble)*dft[k] * exp(-2 * M_PI*i*(double)k*(double)n / numSamplesDouble);
		}
	}

	return signal;
}
