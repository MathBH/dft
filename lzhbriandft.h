#include "complex.h"

/*
	A working DFT algorithm taken from https://github.com/lzhbrian/Fast-Fourier-Transform.git
	Used for testing results of this project's FFT algorithms.
	
	File name changed to "lzhbriandft.h" to differentiate between other dfts.
	
	See below for author information.
	For copyrigth and license, see file "lzhbrian-Fast-Fourier-Transform-LICENSE" in dft folder.
*/

// Work by Lin, Tzu-Heng
// W42, 2014011054
// Dept. of Electronic Engineering, Tsinghua University
// DSP Course Work


complex* DFT(complex input_seq[], int N);

/********************************************************************/
// DFT
	// input_seq[]: 
	// N: size of input_seq
complex* DFT(complex input_seq[], int N) {

	// Calc WN
	complex* WN = new complex[N];
	WN = Calc_WN(N);

	complex* DFTed_seq = new complex[N];
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			int k_mod = (i*j) % N;
			DFTed_seq[i] = ComplexAdd( DFTed_seq[i], ComplexMul(input_seq[j], WN[k_mod]) );
		}
	}

	return DFTed_seq;
}