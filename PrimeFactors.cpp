#include "stdafx.h"
#include "PrimeFactors.h"
#include <math.h>
#include <algorithm>

PrimeFactors::PrimeFactors()
{
}


PrimeFactors::~PrimeFactors()
{
}

/*
A function that returns all prime factors of a given number n
Original code from: https://www.geeksforgeeks.org/print-all-prime-factors-of-a-given-number/
Modified to return a pair of vectors, opne with the prime factors and the next with their corresponding prime bases
*/
std::vector<PrimeFactor> PrimeFactors::primeFactors(int n)
{
	//std::vector<PrimeFactor> primesBuffer = std::vector<PrimeFactor>();

	std::vector<PrimeFactor> result = std::vector<PrimeFactor>();
	result.resize(1, PrimeFactor(1));
	//primesBuffer.resize(1, PrimeFactor(1, 1));
	int lastidx = 0;
	PrimeFactor pfBuffer;
	int baseBuffer;

	// (edit)add 2s to result(edit) that divide n 
	while (n % 2 == 0)
	{
		if (result[lastidx].getBase() == 2)
		{
			result[lastidx].incrementPower();
		}
		else
		{
			result.push_back(PrimeFactor(2));
			lastidx++;
		}
		n = n / 2;
	}

	// n must be odd at this point. So we can skip 
	// one element (Note i = i +2) 
	for (int i = 3; i <= sqrt(n); i = i + 2)
	{
		// While i divides n, (edit)add i to result(edit) and divide n 
		while (n%i == 0)
		{
			if (result[lastidx].getBase() == i)
			{
				result[lastidx].incrementPower();
			}
			else
			{
				result.push_back(PrimeFactor(i));
				lastidx++;
			}
			n = n / i;
		}
	}

	// This condition is to handle the case when n 
	// is a prime number greater than 2 
	if (n > 2)
	{
		result.push_back(PrimeFactor(n));
	}

	std::sort(result.begin(), result.end());
	
	return result;
}