#include "stdafx.h"
#include "PfftKMap.h"


PfftKMap::PfftKMap()
{
}

PfftKMap::PfftKMap(int size, std::vector<PrimeFactor> primeFactors)
{
	N = size;
	Ni = primeFactors;
	d = Ni.size();
	lcSolver = LCSolver();
}


PfftKMap::~PfftKMap()
{
}

/*
	"Fast Fourier Transform and Convolution Algorithms" p.126

	edits: i starts at 0 instead and stops at d-1
*/
int PfftKMap::getIndex(std::vector<int> indicies)
{
	int n = 0;
	int nn;
	int t;
	for (int i = 0; i < d; i++) {
		nn = N / Ni[i].getValue();
		t = lcSolver.solvePfftTi(nn, Ni[i]);
		n += N / Ni[i].getValue()*t*indicies[i];
	}
	n %= N;
	return n;
}
