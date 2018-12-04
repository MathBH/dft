#include "stdafx.h"
#include "PrimeFactors.h"
#include <math.h>
#include <algorithm>

class PrimeFactor
{
public:

	int value;
	int base;

	PrimeFactor(int v, int b)
	{
		value = v;
		base = b;
	}
	PrimeFactor()
	{
		value = 0;
		base = 0;
	}
	~PrimeFactor(){}

	bool operator < (const PrimeFactor& pf) const
	{
		return value < pf.value;
	}
};

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
std::vector<int> PrimeFactors::primeFactors(int n)
{
	//std::vector<PrimeFactor> primesBuffer = std::vector<PrimeFactor>();

	std::vector<int> result = std::vector<int>();
	result.resize(1, 1);
	//primesBuffer.resize(1, PrimeFactor(1, 1));
	int lastidx = 0;

	// (edit)add 2s to result(edit) that divide n 
	while (n % 2 == 0)
	{
		if (result[lastidx]%2 == 0)
		{
			result[lastidx] *= 2;
		}
		else
		{
			result.push_back(2);
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
			if (result[lastidx]%i== 0)
			{
				result[lastidx] *= i;
			}
			else
			{
				result.push_back(i);
				lastidx++;
			}
			n = n / i;
		}
	}

	// This condition is to handle the case when n 
	// is a prime number greater than 2 
	if (n > 2)
	{
		result.push_back(n);
	}

	std::sort(result.begin(), result.end());
	
	return result;
}