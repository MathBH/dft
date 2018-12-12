#include "stdafx.h"
#include "CooleyTukeyFFT.h"
#define _USE_MATH_DEFINES
#include <math.h>


CooleyTukeyFFT::CooleyTukeyFFT()
{
}


CooleyTukeyFFT::~CooleyTukeyFFT()
{
}

std::vector<std::complex<double>> CooleyTukeyFFT::getDft(std::vector<std::complex<double>> samples)
{
	int n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;

	// change samples to book format - TODO: refactor this omg pls
	int samplesCopyLength = samples.size() * 2;
	int tempIdx;
	std::vector<double> samplesCopy = std::vector<double>(samplesCopyLength);
	std::complex<double> tempc;
	for (int idx = 1; idx < samplesCopyLength; idx += 2)
	{
		tempIdx = idx * 0.5;
		tempc = samples[tempIdx];
		samplesCopy[idx] = tempc.imag();
		samplesCopy[idx - 1] = tempc.real();
	}
	// change samples to book format

	int nn = samplesCopyLength /2;
	n = nn << 1;
	j = 1;

	// Bit-Reversal
	for (i = 1; i < n; i += 2)
	{
		if (j > i)
		{
			std::swap(samplesCopy[j - 1], samplesCopy[i - 1]);
			std::swap(samplesCopy[j], samplesCopy[i]);
		}
		m = nn;
		while (m >= 2 && j > m)
		{
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	// Danielson-Lanczos
	mmax = 2;
	while (n > mmax)
	{
		istep = mmax << 1;
		theta = (2.*M_PI) / (double)mmax;
		wtemp = sin(0.5*theta);
		wpr = -2.*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.;
		wi = 0.;
		for (m = 1; m < mmax; m += 2)
		{
			for (i = m; i <= n; i += istep)
			{
				j = i + mmax;
				tempr = wr * samplesCopy[j - 1] - wi * samplesCopy[j];
				tempi = wr * samplesCopy[j] + wi * samplesCopy[j - 1];
				samplesCopy[j - 1] = samplesCopy[i - 1] - tempr;
				samplesCopy[j] = samplesCopy[i] - tempi;
				samplesCopy[i - 1] += tempr;
				samplesCopy[i] += tempi;
			}
			wr = (wtemp = wr)*wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}

	// put in complex<double> vector - TODO: refactor this omg pls
	std::complex<double> imagVal = -1;
	imagVal = sqrt(imagVal);
	std::vector<std::complex<double>> dft = std::vector<std::complex<double>>(samplesCopyLength / 2);
	for (int idx = 1; idx < samplesCopyLength; idx += 2)
	{
		tempIdx = idx * 0.5;
		dft[tempIdx] = -imagVal*samplesCopy[idx] + samplesCopy[idx - 1];
	}
	// put in complex<double> vector

	return dft;
}
