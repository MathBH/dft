#pragma once
#include "DFTAlgorithm.h"

/*
	"Fast Fourier Transform and Convolution" p.145 - p150
*/

class DFT2 : public DFTAlgorithm
{
public:
	DFT2() {}
	~DFT2() {}
	std::vector<std::complex<double>> getDft(std::vector<std::complex<double>> samples) override; //TODO: refactor with object for complex list structure
};

class DFT3 : public DFTAlgorithm
{
public:
	DFT3() {}
	~DFT3() {}
	std::vector<std::complex<double>> getDft(std::vector<std::complex<double>> samples) override; //TODO: refactor with object for complex list structure
};

class DFT4 : public DFTAlgorithm
{
public:
	DFT4() {}
	~DFT4() {}
	std::vector<std::complex<double>> getDft(std::vector<std::complex<double>> samples) override; //TODO: refactor with object for complex list structure
};

class DFT5 : public DFTAlgorithm
{
public:
	DFT5() {}
	~DFT5() {}
	std::vector<std::complex<double>> getDft(std::vector<std::complex<double>> samples) override; //TODO: refactor with object for complex list structure
};

class DFT7 : public DFTAlgorithm
{
public:
	DFT7() {}
	~DFT7() {}
	std::vector<std::complex<double>> getDft(std::vector<std::complex<double>> samples) override; //TODO: refactor with object for complex list structure
};