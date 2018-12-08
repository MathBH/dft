#pragma once
#include <vector>

class PrimeFactor
{
private:
	int value;
	int base;
	int power;

public:

	PrimeFactor(int b, int p)
	{
		base = b;
		power = p;
		value = (int) pow(base, power);
	}

	PrimeFactor(int b)
	{
		base = b;
		power = 1;
		value = base;
	}

	PrimeFactor()
	{
		value = 0;
		base = 0;
	}

	~PrimeFactor() {}

	int getValue() { return value; }
	int getBase() { return base; }
	int getPower() { return power; }

	void incrementPower()
	{
		power++;
		value *= base;
	}

	bool operator < (const PrimeFactor& pf) const
	{
		return value < pf.value;
	}
};

class PrimeFactors
{
public:
	PrimeFactors();
	~PrimeFactors();
	std::vector<PrimeFactor> primeFactors(int n);
};

