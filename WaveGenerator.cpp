#include "stdafx.h"
#include "WaveGenerator.h"
#define _USE_MATH_DEFINES
#include <math.h>

#define TWO_PI (2.0*M_PI)

WaveGenerator::WaveGenerator()
{
}


WaveGenerator::~WaveGenerator()
{
}

std::vector<std::complex<double>> WaveGenerator::sinwave(double frequency, int numSamples, double samplingPeriod, double amplitude)
{
	double increment = TWO_PI*frequency*samplingPeriod;

	std::vector<std::complex<double>> sinwave = std::vector<std::complex<double>>(numSamples);
	for (int i = 0; i < numSamples; i++) {
		sinwave[i] = amplitude*sin((double)i*increment);
	}

	return sinwave;
}

std::vector<std::complex<double>> WaveGenerator::coswave(double frequency, int numSamples, double samplingPeriod, double amplitude)
{
	double increment = TWO_PI*frequency*samplingPeriod;

	std::vector<std::complex<double>> coswave = std::vector<std::complex<double>>(numSamples);
	for (int i = 0; i < numSamples; i++) {
		coswave[i] = amplitude*cos((double)i*increment);
	}

	return coswave;
}
