#pragma once
#include <vector>
#include <complex>

class WaveGenerator
{
public:
	WaveGenerator();
	~WaveGenerator();
	std::vector<std::complex<double>> sinwave(double frequency, int numSamples, double samplingPeriod, double amplitude);
	std::vector<std::complex<double>> coswave(double frequency, int numSamples, double samplingPeriod, double amplitude);
};

