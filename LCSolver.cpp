#include "stdafx.h"
#include "LCSolver.h"
#include <cmath>

LCSolver::LCSolver()
{
}


LCSolver::~LCSolver()
{
}

/*
	Mod operation for handling large powers.
	Taken from https://stackoverflow.com/questions/30244306/find-the-modulo-of-division-of-very-big-numbers
*/
int largeMod(int pow, int val, int MOD)
{
	if (pow == 0)
		return 1;
	int v = largeMod(pow / 2, val, MOD);
	if (pow % 2 == 0)
		return (v*v) % MOD;
	else
		return (((v*val) % MOD) * v) % MOD;
}

double totient(PrimeFactor pf)
{
	double value = pf.getValue();
	if (pf.getPower() == 1)
		return value - 1.;

	double base = pf.getBase();
	return value * (1. - (1. / base));
}

double LCSolver::solvePfftTi(int a, PrimeFactor m)
{
	int pow = (int)totient(m) - 1;
	int result = largeMod(pow, a, m.getValue());
	return result;
}