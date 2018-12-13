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
	std::complex<double> _i = -1.;
	_i = sqrt(_i);

	std::vector<std::complex<double>> ctdft = std::vector<std::complex<double>>(samples);

	int n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;
	n = ctdft.size();
	j = 0;

	// Bit-Reversal				// All in all, fairly familiar to the prime-factor splitting into multiples
	for (i = 0; i <= n; i++)	// bassically placing values in output order from odd-index/evem-index split each time
	{
		if (j > i)
		{
			std::swap(ctdft[j], ctdft[i]);
		}
		m = n >> 1;
		while (m >= 1 && j >= m)
		{
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	// Danielson-Lanczos
	mmax = 1;
	while (n > mmax)
	{
		// every loop of while is multiplying mmax by 2, effectively switching dimensions each time (if we see this as a multidimensional thing with each dimension of size 2)
		istep = mmax << 1;
		theta = M_PI / (double)mmax;		// These are the values for at a given dimension (because row*col becomes new row length so then vals dependent on N use mmax)
		wtemp = sin(0.5*theta);					// notice theta uses 2*pi, coz essentialyl a 2-point dft
		wpr = -2.*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.;
		wi = 0.;
		for (m = 0; m < mmax; m++) // is bassically iterating start pos through the increment like with when u were mapping Prime-Factor
		{							  // start pos increases by 2 because of the imaginary/real representation
			for (i = m; i < n; i += istep) // So istep (the increment) starts at 4, then 8, then 16, 32, etc (ie: 2 exponential)
			{								// So, bassically iteratino through the whole thing at increasing power 2 intervals
				// And this is the actual Danielson-Lancoz calculation part
				j = i + mmax;
				tempr = wr * ctdft[j].real() - wi * -ctdft[j].imag();	// 4 entries coz image and real
				tempi = wr * -ctdft[j].imag() + wi * ctdft[j].real();

				ctdft[j] = ctdft[i] - (tempr - _i*tempi);
				ctdft[i] += (tempr - _i*tempi);
			}
			wtemp = wr;
			wr = wtemp*wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}

	return ctdft;
}
