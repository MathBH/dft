#include "stdafx.h"
#include "PFFT.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include "PfftKMap.h"
#include "PfftNMap.h"
#include "MultiDimensionalCounter.h"

/*
	Implementation of PFFT algorithm from:
	"Fast Fourier Transform and Convolution Algorithms" p.125-129
*/

PFFT::PFFT()
{
	pfCalulator = PrimeFactors();
}


PFFT::~PFFT()
{
}

std::vector<std::complex<double>> PFFT::getDft(std::vector<std::complex<double>> samples)
{
	std::complex<double> i = -1; // define imaginary number
	i = sqrt(i);

	int numSamples = samples.size();	// calculate numSamples and its prime factors
	std::vector<PrimeFactor> primeFactors = pfCalulator.primeFactors(numSamples);

	//PfftNMap nMap = PfftNMap(numSamples, primeFactors); // build multidimensional n mapping
	//PfftKMap kMap = PfftKMap(numSamples, primeFactors); // build multidimensional k mapping
	//
	//int numDimensions = primeFactors.size();
	//std::vector<int> dimensions = std::vector<int>(numDimensions);
	//for (int i = 0; i < numDimensions; i++)
	//{
	//	dimensions[numDimensions - 1 - i] = primeFactors[i].getValue();
	//}
	//if (dimensions.size() < 1) {
	//	dimensions.push_back(1);
	//}
	//MultiDimensionalCounter mcd = MultiDimensionalCounter(dimensions);

	std::vector<std::complex<double>> dft = std::vector<std::complex<double>>(numSamples);

	//std::vector<int> idx;
	//int n;
	//int k;
	//while (mcd.hasNext()) {
	//	idx = mcd.nextPos();
	//	n = nMap.getIndex(idx);
	//	k = kMap.getIndex(idx);
	//	dft[k] += samples[n] * exp(2 * M_PI*-i*(double)k*(double)n / (double)numSamples);
	//}
	//Todo: edit this to iterate through multidimensional indicies instead and do the ops
	//for (int n = 0; n < numSamples; n++)
	//{
	//	for (int k = 0; k < numSamples; k++)
	//	{
	//		dft[n] += samples[k] * exp(2 * M_PI*i*(double)k*(double)n / (double)numSamples);
	//	}
	//}

	return dft;
}
