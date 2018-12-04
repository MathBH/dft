#pragma once
#include <vector>

class WaveGenerator
{
public:
	WaveGenerator();
	~WaveGenerator();
	std::vector<double> sinwave(double frequency, int numSamples, double samplingPeriod, double amplitude);
	std::vector<double> coswave(double frequency, int numSamples, double samplingPeriod, double amplitude);
};

