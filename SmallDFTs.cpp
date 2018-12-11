#include "stdafx.h"
#include "SmallDFTs.h"
#define _USE_MATH_DEFINES
#include <math.h>

/*
	2-Point Dft
	"Fast Fourier Transform and Convolution" p.145
*/
std::vector<std::complex<double>> DFT2::getDft(std::vector<std::complex<double>> samples)
{
	std::vector<std::complex<double>> X = std::vector<std::complex<double>>(2);

	std::complex<double> m0 = samples[0] + samples[1];
	std::complex<double> m1 = samples[0] - samples[1];
	
	X[0] = m0;
	X[1] = m1;

	return X;
}

/*
	3-Point Dft
	"Fast Fourier Transform and Convolution" p.145
*/
std::vector<std::complex<double>> DFT3::getDft(std::vector<std::complex<double>> samples)
{
	std::complex<double> j = -1.; //TODO: refactor coz it's pasted everywhere
	j = sqrt(j);

	std::vector<std::complex<double>> X = std::vector<std::complex<double>>(3);

	std::complex<double> u = 2. * M_PI/3.;
	std::complex<double> t1 = samples[1] + samples[2];
	std::complex<double> m0 = samples[0] + t1;
	std::complex<double> m1 = (cos(u) - 1.)*t1;
	std::complex<double> m2 = j * sin(u)*(samples[2] - samples[1]);
	std::complex<double> s1 = m0 + m1;

	X[0] = m0;
	X[1] = s1 + m2;
	X[2] = s1 - m2;

	return X;
}

/*
4-Point Dft
"Fast Fourier Transform and Convolution" p.145-146
*/
std::vector<std::complex<double>> DFT4::getDft(std::vector<std::complex<double>> samples)
{
	std::complex<double> j = -1.; //TODO: refactor coz it's pasted everywhere
	j = sqrt(j);

	std::vector<std::complex<double>> X = std::vector<std::complex<double>>(4);

	std::complex<double> t1 = samples[0] + samples[2];
	std::complex<double> t2 = samples[1] + samples[3];
	std::complex<double> m0 = t1 + t2;
	std::complex<double> m1 = t1 - t2;
	std::complex<double> m2 = samples[0] - samples[2];
	std::complex<double> m3 = j*(samples[3] - samples[1]);

	X[0] = m0;
	X[1] = m2 + m3;
	X[2] = m1;
	X[3] = m2 - m3;

	return X;
}

/*
5-Point Dft
"Fast Fourier Transform and Convolution" p.146
*/
std::vector<std::complex<double>> DFT5::getDft(std::vector<std::complex<double>> samples)
{
	std::complex<double> j = -1.; //TODO: refactor coz it's pasted everywhere
	j = sqrt(j);

	std::vector<std::complex<double>> X = std::vector<std::complex<double>>(5);

	std::complex<double> u = 2. * M_PI / 5.;
	std::complex<double> t1 = samples[1] + samples[4];
	std::complex<double> t2 = samples[2] + samples[3];
	std::complex<double> t3 = samples[1] - samples[4];
	std::complex<double> t4 = samples[3] - samples[2];
	std::complex<double> t5 = t1 + t2;
	std::complex<double> m0 = 1.*(samples[0] + t5);
	std::complex<double> m1 = ((cos(u) + cos(2.*u)) / 2. - 1.)*t5;
	std::complex<double> m2 = ((cos(u) - cos(2.*u)) / 2.)*(t1 - t2);
	std::complex<double> m3 = -j * (sin(u)) * (t3 + t4);
	std::complex<double> m4 = -j * (sin(u) + sin(2.*u)) * t4;
	std::complex<double> m5 = j * (sin(u) - sin(2.*u))*t3;
	std::complex<double> s3 = m3 - m4;
	std::complex<double> s5 = m3 + m5;
	std::complex<double> s1 = m0 + m1;
	std::complex<double> s2 = s1 + m2;
	std::complex<double> s4 = s1 - m2;

	X[0] = m0;
	X[1] = s2 + s3;
	X[2] = s4 + s5;
	X[3] = s4 - s5;
	X[4] = s2 - s3;

	return X;
}

/*
7-Point Dft
"Fast Fourier Transform and Convolution" p.146-147
*/
std::vector<std::complex<double>> DFT7::getDft(std::vector<std::complex<double>> samples)
{
	std::complex<double> j = -1; //TODO: refactor coz it's pasted everywhere
	j = sqrt(j);

	std::vector<std::complex<double>> X = std::vector<std::complex<double>>(7);

	std::complex<double> u = 2. * M_PI / 7.;

	std::complex<double> t1 = samples[1] + samples[6];
	std::complex<double> t2 = samples[2] + samples[5];
	std::complex<double> t3 = samples[3] + samples[4];
	std::complex<double> t4 = t1 + t2 + t3;
	std::complex<double> t5 = samples[1] - samples[6];
	std::complex<double> t6 = samples[2] - samples[5];
	std::complex<double> t7 = samples[4] - samples[3];
	std::complex<double> t8 = t1 - t3;
	std::complex<double> t9 = t3 - t2;
	std::complex<double> t10 = t5 + t6 + t7;
	std::complex<double> t11 = t7 - t5;
	std::complex<double> t12 = t6 - t7;

	std::complex<double> m0 = samples[0] + t4;
	std::complex<double> m1 = ((cos(u) + cos(2.*u) + cos(3.*u)) / 3. - 1.)*t4;

	std::complex<double> t13 = -t8 - t9;

	std::complex<double> m2 = ((2.*cos(u) - cos(2.*u) - cos(3.*u))/3.) *t8 ;
	std::complex<double> m3 = (( cos(u) - 2.*cos(2.*u) + cos(3.*u))/3.)*t9 ;
	std::complex<double> m4 = (( cos(u) + cos(2.*u) - 2.*cos(3.*u))/3.)*t13;

	std::complex<double> s0 = -m2 - m3;
	std::complex<double> s1 = -m2 - m4;

	std::complex<double> m5 = -j * ((sin(u) + sin(2.*u) - sin(3.*u)) / 3.)*t10;

	std::complex<double> t14 = -t11 - t12;

	std::complex<double> m6 = j*((2.*sin(u) - sin(2.*u) + sin(3.*u))/3.)*t11;
	std::complex<double> m7 = j*((sin(u) - 2.*sin(2.*u) - sin(3.*u))/3.)*t12;
	std::complex<double> m8 = j*((sin(u) + sin(2.*u) + 2.*sin(3.*u))/3.)*t14;

	std::complex<double> s2 = -m6 - m7;
	std::complex<double> s3 = m6 + m8;
	std::complex<double> s4 = m0 + m1;
	std::complex<double> s5 = s4 - s0;
	std::complex<double> s6 = s4 + s1;
	std::complex<double> s7 = s4 + s0 - s1;
	std::complex<double> s8 = m5 - s2;
	std::complex<double> s9 = m5 - s3;
	std::complex<double> s10 = m5 + s2 + s3;

	X[0] = m0;
	X[1] = s5 + s8;
	X[2] = s6 + s9;
	X[3] = s7 - s10;
	X[4] = s7 + s10;
	X[5] = s6 - s9;
	X[6] = s5 - s8;

	return X;
}
