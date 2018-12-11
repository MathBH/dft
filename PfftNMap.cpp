#include "stdafx.h"
#include "PfftNMap.h"


PfftNMap::PfftNMap()
{
}

PfftNMap::PfftNMap(int size, std::vector<PrimeFactor> primeFactors)
{
	N = size;
	Ni = primeFactors;
	d = Ni.size();
}


PfftNMap::~PfftNMap()
{
}

/*
	"Fast Fourier Transform and Convolution Algorithms" p.126

	edits: i starts at 0 instead and stops at d-1
*/
int PfftNMap::getIndex(std::vector<int> indicies)
{
	int n = 0;
	for (int i = 0; i < d; i++) {
		n += N / Ni[i].getValue()*indicies[i];
	}
	n %= N;
	return n;
}
