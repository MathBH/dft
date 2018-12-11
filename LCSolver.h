#pragma once
#include "PrimeFactors.h"
/*
	For solving linear congruences
	
	LCSolver is currently only adapted for solving for ti in PFFT actually so maybe rename..
*/

class LCSolver
{
protected:
public:
	LCSolver();
	~LCSolver();
	double solvePfftTi(int a, PrimeFactor m);
};

