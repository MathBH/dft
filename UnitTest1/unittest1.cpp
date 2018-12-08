#include "stdafx.h"
#include "CppUnitTest.h"
#include "..\PFFTMapper.h"
#include "..\PFFTMapper.cpp"
#include "..\PrimeFactors.h"
#include "..\PrimeFactors.cpp"
#include "..\WaveGenerator.h"
#include "..\WaveGenerator.cpp"
#include "..\PrimeFactorDft.h"
#include "..\PrimeFactorDft.cpp"
#include "..\BruteForceDft.h"
#include "..\BruteForceDft.cpp"
#include "..\SmallDfts.h"
#include "..\SmallDfts.cpp"
#include "..\DFTAlgorithm.h"
#include <vector>
#include <iostream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest1
{		
	TEST_CLASS(PFFTMapperTests)
	{
	public:
		PFFTMapper pfftMapper = PFFTMapper();
		TEST_METHOD(MappingTest1)
		{
			int numValues = 12;
			std::vector<PrimeFactor> primes = std::vector<PrimeFactor>({ PrimeFactor(1),PrimeFactor(3),PrimeFactor(2,2) });
			std::vector<int> valuesMapped = pfftMapper.basicMapping(numValues, primes);
			std::vector<int> expectedMapping = std::vector<int>({ 0,4,8,3,7,11,6,10,2,9,1,5 });
			Assert::IsTrue(valuesMapped.size() == expectedMapping.size());
			for (int i = 0; i < numValues; i++)
			{
				Assert::IsTrue(valuesMapped[i] == expectedMapping[i]);
			}
		}
		TEST_METHOD(MappingTest2)
		{
			int numValues = 30;
			std::vector<PrimeFactor> primes = std::vector<PrimeFactor>({ PrimeFactor(1),PrimeFactor(2),PrimeFactor(3),PrimeFactor(5) });
			std::vector<int> valuesMapped = pfftMapper.basicMapping(numValues, primes);
			std::vector<int> expectedMapping = std::vector<int>({ 0,15,10,25,20,5,6,21,16,1,26,11,12,27,22,7,2,17,18,3,28,13,8,23,24,9,4,19,14,29 });
			Assert::IsTrue(valuesMapped.size() == expectedMapping.size());
			for (int i = 0; i < numValues; i++)
			{
				Assert::IsTrue(valuesMapped[i] == expectedMapping[i]);
			}
		}
	};

	TEST_CLASS(PrimeFactorsTests)
	{
	private:
		void assertSameNumbers(std::vector<PrimeFactor> A, std::vector<PrimeFactor> B)
		{
			int aLen = A.size();
			int bLen = B.size();
			Assert::IsTrue(aLen == bLen);
			for (int i = 0; i < aLen; i++)
			{
				Assert::IsTrue(A[i].getValue() == B[i].getValue());
				Assert::IsTrue(A[i].getBase() == B[i].getBase());
				Assert::IsTrue(A[i].getPower() == B[i].getPower());
			}
		}
	public:
		int val;
		std::vector<PrimeFactor> result;
		std::vector<PrimeFactor> primesExpected;
		PrimeFactors pf = PrimeFactors();
		TEST_METHOD(PrimeOf12)
		{
			val = 12;
			primesExpected = std::vector<PrimeFactor>({ PrimeFactor(1),PrimeFactor(3),PrimeFactor(2,2) });
			result = pf.primeFactors(val);
			assertSameNumbers(result, primesExpected);

		}
		TEST_METHOD(PrimeOf30)
		{
			val = 30;
			primesExpected = std::vector<PrimeFactor>({ PrimeFactor(1),PrimeFactor(2),PrimeFactor(3),PrimeFactor(5) });
			result = pf.primeFactors(val);
			assertSameNumbers(result, primesExpected);
		}
	};

	TEST_CLASS(PFFTMapperAndPrimeFactorsTests)
	{
	private:
		void assertSameNumbers(std::vector<int> A, std::vector<int> B)
		{
			int aLen = A.size();
			int bLen = B.size();
			Assert::IsTrue(aLen == bLen);
			for (int i = 0; i < aLen; i++)
			{
				Assert::IsTrue(A[i] == B[i]);
			}
		}
	public:
		PFFTMapper pfftMapper = PFFTMapper();
		PrimeFactors pf = PrimeFactors();
		TEST_METHOD(MappingTest1)
		{
			int numValues = 12;
			std::vector<PrimeFactor> primes = pf.primeFactors(numValues);
			std::vector<int> valuesMapped = pfftMapper.basicMapping(numValues, primes);
			std::vector<int> expectedMapping = std::vector<int>({ 0,4,8,3,7,11,6,10,2,9,1,5 });
			Assert::IsTrue(valuesMapped.size() == expectedMapping.size());
			for (int i = 0; i < numValues; i++)
			{
				Assert::IsTrue(valuesMapped[i] == expectedMapping[i]);
			}
		}
		TEST_METHOD(MappingTest2)
		{
			int numValues = 30;
			std::vector<PrimeFactor> primes = pf.primeFactors(numValues);
			std::vector<int> valuesMapped = pfftMapper.basicMapping(numValues, primes);
			std::vector<int> expectedMapping = std::vector<int>({ 0,15,10,25,20,5,6,21,16,1,26,11,12,27,22,7,2,17,18,3,28,13,8,23,24,9,4,19,14,29 });
			Assert::IsTrue(valuesMapped.size() == expectedMapping.size());
			for (int i = 0; i < numValues; i++)
			{
				Assert::IsTrue(valuesMapped[i] == expectedMapping[i]);
			}
		}
	};


	TEST_CLASS(PFFTTest)
	{
	private:
		void assertSimilar(std::vector<std::complex<double>> A, std::vector<std::complex<double>> B, double errorThreshold)
		{
			int aLen = A.size();
			int bLen = B.size();
			double error;
			Assert::IsTrue(aLen == bLen);
			for (int i = 0; i < aLen; i++)
			{
				Assert::IsTrue(abs(A[i].imag() - B[i].imag()) < errorThreshold);
				Assert::IsTrue(abs(A[i].real() - B[i].real()) < errorThreshold);
			}
		}
	public:
		WaveGenerator waveGen = WaveGenerator();
		PrimeFactorDft pfft = PrimeFactorDft();
		BruteForceDft bdft = BruteForceDft();
		std::vector<std::complex<double>> dftA;
		std::vector<std::complex<double>> dftB;
		std::vector<std::complex<double>> wave;
		double error;
		TEST_METHOD(CosWaveTest30)
		{
			error = 0.000000001;
			wave = waveGen.coswave(10., 30, 0.01, 1.);
			dftA = bdft.getDft(wave);
			dftB = pfft.getDft(wave);
			assertSimilar(dftA, dftB, error);
		}
		TEST_METHOD(SinWaveTest30)
		{
			error = 0.000000001;
			wave = waveGen.sinwave(10., 30, 0.01, 1.);
			dftA = bdft.getDft(wave);
			dftB = pfft.getDft(wave);
			assertSimilar(dftA, dftB, error);
		}
		TEST_METHOD(CosWaveTest12)
		{
			error = 0.000000001;
			wave = waveGen.coswave(10., 12, 0.01, 1.);
			dftA = bdft.getDft(wave);
			dftB = pfft.getDft(wave);
			assertSimilar(dftA, dftB, error);
		}
		TEST_METHOD(SinWaveTest12)
		{
			error = 0.000000001;
			wave = waveGen.sinwave(10., 12, 0.01, 1.);
			dftA = bdft.getDft(wave);
			dftB = pfft.getDft(wave);
			assertSimilar(dftA, dftB, error);
		}
	};
}