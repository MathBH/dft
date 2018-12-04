#include "stdafx.h"
#include "PrimeFactorDft.h"

PrimeFactorDft::PrimeFactorDft()
{
	mapper = PFFTMapper();
	primeFactorCalculator = PrimeFactors();
}


PrimeFactorDft::~PrimeFactorDft()
{
}

std::vector<std::complex<double>> PrimeFactorDft::getDft(std::vector<double> samples)
{
	int numSamples = samples.size();
	std::vector<int> primeFactors = this->primeFactorCalculator.primeFactors(numSamples);

	int numPrimes = primeFactors.size();
	std::vector<std::vector<std::complex<double>>> db = std::vector<std::vector<std::complex<double>>>(numPrimes,std::vector<std::complex<double>>(numSamples)); // TODO: probably make into multidimentional matrix and move through

	std::vector<int> convolutionArray = this->mapper.basicMapping(numSamples, primeFactors);

	for (int i = 0; i < numSamples; i++)
	{
		db[0][i] = samples[convolutionArray[i]];
	}

	// Compute DFT
	// Notice each row of db is effectively a multidimensional matrix mapping.
	// We iterate through it by assigning to numVectors the number of vectors in dimension d,
	// to increment the value that takes us to the element in a given vector, and to vectorLength
	// the length of the vectors in dimension d.
	// 
	// We then compute the dfts for the vectors in dimension d, shift to the next dimension, and
	// proceed until all dimensions have been navigated.
	//
	// Vectors to perform the fourier transforms on are defined by assigning an index value to
	// vectorStart and vectorEnd representing the first and last element of the vector.
	//
	// We use a hashtable to select the right dft algorithm for the vector.
	//
	// We then iterate over the defined vector and compute the dft into db[d]

	int j;
	int vectorLength;
	int vectorOffset;
	int increment = 1;

	int vectorStart;
	int vectorEnd;

	int currPos;
	std::vector<std::complex<double>> vectorBuffer;
	for (int d = 1; d < numPrimes; d++)
	{
		increment *= primeFactors[d - 1];
		vectorLength = primeFactors[d];
		vectorOffset = vectorLength * increment;

		// use primeFactors[d] to select dft alg
		vectorBuffer = std::vector<std::complex<double>>(vectorLength);

		// for each value betwen increments
		for (int j = 0; j < increment; j++)
		{
			// for vectors in that increment
			for (int i = j; i < numSamples; i += vectorOffset)
			{
				// build vector for dft alg to process
				for (int h = 0; h < vectorLength; h++)
				{
					currPos = i + h * increment;
					vectorBuffer[h] = db[d-1][currPos];
				}
				// feed dft vals into current db
				for (int h = 0; h < vectorLength; h++)
				{
					currPos = i + h * increment;
					db[d][currPos] = h; // TODO compute dft for X_h of vectorBuffer
				}
			}
		}
	}

	// TODO: reverse convolution

	return db[numPrimes-1];
}
