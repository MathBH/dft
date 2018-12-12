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
#include "..\PFFT.h"
#include "..\PFFT.cpp"
#include "..\SmallDfts.h"
#include "..\SmallDfts.cpp"
#include "..\DFTAlgorithm.h"
#include "..\LCSolver.h"
#include "..\LCSolver.cpp"
#include "..\PfftKMap.h"
#include "..\PfftKMap.cpp"
#include "..\PfftNMap.h"
#include "..\PfftNMap.cpp"
#include "..\PFFTMap.h"
#include "..\MultiDimensionalCounter.cpp"
#include "..\MultiDimensionalCounter.h"
#include "..\CooleyTukeyFFT.cpp"
#include "..\CooleyTukeyFFT.h"
#include <vector>
#include <iostream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest1
{		
	TEST_CLASS(PFFTMapperTests)
	{
	public:
		PFFTMapper pfftMapper = PFFTMapper();
		TEST_METHOD(InputMappingTest1)
		{
			int numValues = 12;
			std::vector<PrimeFactor> primes = std::vector<PrimeFactor>({ PrimeFactor(1),PrimeFactor(3),PrimeFactor(2,2) });
			std::vector<int> valuesMapped = pfftMapper.inputMapping(numValues, primes);
			std::vector<int> expectedMapping = std::vector<int>({ 0,4,8,3,7,11,6,10,2,9,1,5 });
			Assert::IsTrue(valuesMapped.size() == expectedMapping.size());
			for (int i = 0; i < numValues; i++)
			{
				Assert::IsTrue(valuesMapped[i] == expectedMapping[i]);
			}
		}
		TEST_METHOD(InputMappingTest2)
		{
			int numValues = 30;
			std::vector<PrimeFactor> primes = std::vector<PrimeFactor>({ PrimeFactor(1),PrimeFactor(2),PrimeFactor(3),PrimeFactor(5) });
			std::vector<int> valuesMapped = pfftMapper.inputMapping(numValues, primes);
			std::vector<int> expectedMapping = std::vector<int>({ 0,15,10,25,20,5,6,21,16,1,26,11,12,27,22,7,2,17,18,3,28,13,8,23,24,9,4,19,14,29 });
			Assert::IsTrue(valuesMapped.size() == expectedMapping.size());
			for (int i = 0; i < numValues; i++)
			{
				Assert::IsTrue(valuesMapped[i] == expectedMapping[i]);
			}
		}
		TEST_METHOD(OutputMappingTest1)
		{
			int numValues = 12;
			std::vector<PrimeFactor> primes = std::vector<PrimeFactor>({ PrimeFactor(1),PrimeFactor(3),PrimeFactor(2,2) });
			std::vector<int> valuesMapped = pfftMapper.outputMapping(numValues, primes);
			std::vector<int> expectedMapping = std::vector<int>({ 0,4,8,9,1,5,6,10,2,3,7,11 });
			Assert::IsTrue(valuesMapped.size() == expectedMapping.size());
			for (int i = 0; i < numValues; i++)
			{
				Assert::IsTrue(valuesMapped[i] == expectedMapping[i]);
			}
		}
	};

	TEST_CLASS(DFTTest)
	{
	private:
		void assertDoubleEqual(double a, double b, double error)
		{
			Assert::IsTrue(abs(a - b)< error);
		}
	
		// Assume inputs "frequencies" and "amplitudes" have the same length
		void assertDFTMatch(std::vector<std::complex<double>> dft, std::vector<double> cosFrequencies, std::vector<double> cosAmplitudes, std::vector<double> sinFrequencies, std::vector<double> sinAmplitudes, double error, double samplingPeriod)
		{
			std::complex<double> j = -1;
			j = sqrt(j);
	
			std::vector<std::complex<double>> dftCopy = std::vector<std::complex<double>>(dft);
			int numSinFrequencies = sinFrequencies.size();
			int numCosFrequencies = cosFrequencies.size();
			int dftLength = dftCopy.size();
			std::vector<double>::iterator it;
			int f;
			int fidx;
	
			double ampScaler = 1./((double)dftLength / 2.);
			double freqScaler = (double)dftLength * samplingPeriod;

			//double avgAmplitude = 0.;
	
			// get frequencies from reverse second half of dft (positive frequencies)
			for (int i = 0; i < numSinFrequencies; i++) {
				f = sinFrequencies[i];
				fidx = dftLength - f* freqScaler;
				assertDoubleEqual(dftCopy[fidx].imag()*ampScaler, sinAmplitudes[i], error);
				//avgAmplitude += dftCopy[fidx].imag();
				dftCopy[fidx] -= j* dftCopy[fidx].imag(); // zero the imaginary part
			}
			for (int i = 0; i < numCosFrequencies; i++) {
				f = cosFrequencies[i];
				fidx = dftLength - f* freqScaler;
				assertDoubleEqual(dftCopy[fidx].real()*ampScaler, cosAmplitudes[i], error);
				//avgAmplitude += dftCopy[fidx].real();
				dftCopy[fidx] -= dftCopy[fidx].real(); // zero the real part
			}

			//avgAmplitude /= numSinFrequencies;
			//avgAmplitude /= numCosFrequencies;
	
			// iterate backwards to only half of dft (positive frequencies)
			for (int i = dftLength - 1; i > dftLength/2; i--)
			{
				assertDoubleEqual(dftCopy[i].imag(), 0., error);
				assertDoubleEqual(dftCopy[i].real(), 0., error);
			}
		}

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
		// As long as the brute force dft works, the following tests show all other dfts work
		WaveGenerator waveGen = WaveGenerator();
		BruteForceDft bdft = BruteForceDft();
		std::vector<std::complex<double>> wave;
		std::vector<std::complex<double>> rebuiltWave;
		std::vector<std::complex<double>> dft;
		std::vector<double> sinFrequencies;
		std::vector<double> cosFrequencies;
		std::vector<double> sinAmplitudes;
		std::vector<double> cosAmplitudes;
		Logger logger = Logger();
		std::stringstream strStream = std::stringstream();
		std::string stringBuffer;
		double error;
		int numSamples;
		double samplingPeriod;
	
		// Assume same length
		std::vector<std::complex<double>> sumWaves(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b)
		{
			int length = a.size();
			std::vector<std::complex<double>> output = std::vector<std::complex<double>>(length);
	
			for (int i = 0; i < length; i++)
			{
				output[i] = a[i] + b[i];
			}
			
			return output;
		}
	
		// Assume frequency and amplitude lists are equal in size
		std::vector<std::complex<double>> generateWave(int length, double samplingPeriod, std::vector<double> cosFrequencies, std::vector<double> cosAmplitudes, std::vector<double> sinFrequencies, std::vector<double> sinAmplitudes)
		{
			wave = std::vector<std::complex<double>>(length);
			int numSins = sinFrequencies.size();
			int numCos = cosFrequencies.size();
			for (int i = 0; i < numSins; i++) {
				wave = sumWaves(wave,waveGen.sinwave(sinFrequencies[i], length, samplingPeriod, sinAmplitudes[i]));
			}
			for (int i = 0; i < numCos; i++) {
				wave = sumWaves(wave, waveGen.coswave(cosFrequencies[i], length, samplingPeriod, cosAmplitudes[i]));
			}
			return wave;
		}

		/*
			A minor stress test for testing if dfts obtained accurately rebuild the originating waves
		*/
		TEST_METHOD(WaveRebuildTest)
		{
			error = 0.00001;
			double frequency = 20;
			double amplitude = 1.;
			numSamples = 100;

			int freqIterations = 13;
			double freqIncrement = 1.;

			int periodIterations = 17;
			double periodScaleDown = 0.837;

			int ampIterations = 3;
			double ampIncrement = 1.37;

			int nIterations = 11;
			int nIncrement = 3;

			int numTests = freqIterations * periodIterations*ampIterations*nIterations;
			int testCounter = 0;

			for (int f = 0; f < freqIterations; f++)
			{
				cosFrequencies = { frequency };
				sinFrequencies = { frequency };
				samplingPeriod = 1. / (frequency*4.);

				for (int p = 0; p < periodIterations; p++)
				{
					numSamples = (((1. / frequency)) / samplingPeriod)*4.;
					for (int a = 0; a < ampIterations; a++)
					{
						cosAmplitudes = { amplitude };
						sinAmplitudes = { amplitude };

						for (int n = 0; n < nIterations; n++)
						{
							wave = generateWave(numSamples, samplingPeriod, cosFrequencies, cosAmplitudes, sinFrequencies, sinAmplitudes);
							dft = bdft.getDft(wave);
							rebuiltWave = bdft.reverseDft(dft);
							assertSimilar(wave, rebuiltWave, error);
							numSamples += nIncrement;
							testCounter++;
						}
						strStream << "Status: " << testCounter << " of " << numTests << " completed.";
						stringBuffer = strStream.str();
						logger.WriteMessage(stringBuffer.c_str());
						strStream.str("");
						strStream.clear();
						amplitude += ampIncrement;
					}
					samplingPeriod *= periodScaleDown;
				}
				frequency += freqIncrement;
			}
		};
	
		/*
			Note: test still clanky - good enough for extensive testing on the data having the 
			"shape" we want but a fully automated stress test would require some polishing of this
		*/
		TEST_METHOD(WaveDataTest)
		{
			error = 0.00001;
			double frequency = 20;
			double amplitude = 1.;
			numSamples = 100;
	
			//samplingPeriod = 1./(frequency*4.);
	
			int freqIterations = 10;
			double freqIncrement = 1.;
	
			int periodIterations = 7;
			double periodScaleDown = 0.5;
	
			int ampIterations = 3;
			double ampIncrement = 1.37;
	
			int numTests = freqIterations * periodIterations*ampIterations;
			int testCounter = 0;
	
			for (int f = 0; f < freqIterations; f++)
			{
				cosFrequencies = { frequency };
				sinFrequencies = { frequency };
				samplingPeriod = 1. / (frequency*4.);
				//logger.WriteMessage

				for (int p = 0; p < periodIterations; p++)
				{
					numSamples = (((1./frequency))/samplingPeriod)*4.;
					for (int a = 0; a < ampIterations; a++)
					{
						cosAmplitudes = { amplitude };
						sinAmplitudes = { amplitude };

						wave = generateWave(numSamples, samplingPeriod, cosFrequencies, cosAmplitudes, sinFrequencies, sinAmplitudes);
						dft = bdft.getDft(wave);
						assertDFTMatch(dft, cosFrequencies, cosAmplitudes, sinFrequencies, sinAmplitudes, error, samplingPeriod);
						dft.clear();
						wave.clear();

						strStream << "Status: " << testCounter << " of " << numTests << " completed.";
						stringBuffer = strStream.str();
						logger.WriteMessage(stringBuffer.c_str());
						strStream.str("");
						strStream.clear();
						amplitude += ampIncrement;
						testCounter++;
					}
					samplingPeriod *= periodScaleDown;
				}
				frequency += freqIncrement;
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
			primesExpected = std::vector<PrimeFactor>({ PrimeFactor(3),PrimeFactor(2,2) });
			result = pf.primeFactors(val);
			assertSameNumbers(result, primesExpected);

		}
		TEST_METHOD(PrimeOf30)
		{
			val = 30;
			primesExpected = std::vector<PrimeFactor>({ PrimeFactor(2),PrimeFactor(3),PrimeFactor(5) });
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
			std::vector<int> valuesMapped = pfftMapper.inputMapping(numValues, primes);
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
			std::vector<int> valuesMapped = pfftMapper.inputMapping(numValues, primes);
			std::vector<int> expectedMapping = std::vector<int>({ 0,15,10,25,20,5,6,21,16,1,26,11,12,27,22,7,2,17,18,3,28,13,8,23,24,9,4,19,14,29 });
			Assert::IsTrue(valuesMapped.size() == expectedMapping.size());
			for (int i = 0; i < numValues; i++)
			{
				Assert::IsTrue(valuesMapped[i] == expectedMapping[i]);
			}
		}
	};

	TEST_CLASS(LinearCongruenceTest)
	{
	private:
		void assertDoubleEqual(double a, double b, double error)
		{
			Assert::IsTrue(abs(a - b)< error);
		}
	public:
		LCSolver lcSolver = LCSolver();
		double error = 0.0000000001;
		TEST_METHOD(PFFTCasesTest)
		{
			double answer = lcSolver.solvePfftTi(4, PrimeFactor(3));
			assertDoubleEqual(1., answer, error);
			answer = lcSolver.solvePfftTi(3, PrimeFactor(2,2));
			assertDoubleEqual(3., answer, error);
			answer = lcSolver.solvePfftTi(13, PrimeFactor(37));
			assertDoubleEqual(20., answer, error);
			answer = lcSolver.solvePfftTi(8, PrimeFactor(37));
			assertDoubleEqual(14., answer, error);
		}
	};

	TEST_CLASS(MultiDimensionalCounterTest)
	{
	private:
		void assertAreEqual(std::vector<int> a, std::vector<int> b) {
			int alen = a.size();
			int blen = b.size();
			Assert::IsTrue(alen == blen);
			for (int i = 0; i < alen; i++)
			{
				Assert::IsTrue(a[i] == b[i]);
			}
		}
	public:
		MultiDimensionalCounter mdc;
		TEST_METHOD(BasicTest)
		{
			std::vector<int> dimensions = {5,6,9,10,5,6};
			mdc = MultiDimensionalCounter(dimensions);
			for (int a = 0; a < dimensions[5]; a++)
			{
				for (int b = 0; b < dimensions[4]; b++)
				{
					for (int c = 0; c < dimensions[3]; c++)
					{
						for (int d = 0; d < dimensions[2]; d++)
						{
							for (int e = 0; e < dimensions[1]; e++)
							{
								for (int f = 0; f < dimensions[0]; f++)
								{
									Assert::IsTrue(mdc.hasNext());
									assertAreEqual({ f,e,d,c,b,a }, mdc.nextPos());
								}
							}

						}
					}

				}
			}
			Assert::IsFalse(mdc.hasNext());
		}
	};
	TEST_CLASS(SmallDftTest)
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
		BruteForceDft bdft = BruteForceDft();
		std::vector<std::complex<double>> dftA;
		std::vector<std::complex<double>> dftB;
		std::vector<std::complex<double>> wave;
		double error;
		int numSamples;
		TEST_METHOD(DFT2TEST)
		{
			DFT2 dft2 = DFT2();
			numSamples = 2;
			error = 0.000000001;
			int base = 1;
			int cap = 70;
			int increment = 1;
			double a = 10.;
			double c = 0.01;
			double d = 1.;
			for (int N = base; N < cap; N += increment)
			{
				wave = waveGen.sinwave(a, numSamples, c, d);
				dftA = bdft.getDft(wave);
				dftB = dft2.getDft(wave);
				assertSimilar(dftA, dftB, error);

				wave = waveGen.coswave(a, numSamples, c, d);
				dftA = bdft.getDft(wave);
				dftB = dft2.getDft(wave);
				assertSimilar(dftA, dftB, error);

				a += 7.;
				c += 0.009;
				d += 0.05;
			}
		}
		TEST_METHOD(DFT3TEST)
		{
			DFT3 dft3 = DFT3();
			numSamples = 3;
			error = 0.000000001;
			int base = 1;
			int cap = 100;
			int increment = 1;
			double a = 10.;
			double c = 0.01;
			double d = 1.;
			for (int N = base; N < cap; N += increment)
			{
				wave = waveGen.sinwave(a, numSamples, c, d);
				dftA = bdft.getDft(wave);
				dftB = dft3.getDft(wave);
				assertSimilar(dftA, dftB, error);

				wave = waveGen.coswave(a, numSamples, c, d);
				dftA = bdft.getDft(wave);
				dftB = dft3.getDft(wave);
				assertSimilar(dftA, dftB, error);

				a += 7.;
				c += 0.009;
				d += 0.05;
			}
		}
		TEST_METHOD(DFT4TEST)
		{
			DFT4 dft4 = DFT4();
			numSamples = 4;
			error = 0.000000001;
			int base = 1;
			int cap = 100;
			int increment = 1;
			double a = 10.;
			double c = 0.01;
			double d = 1.;
			for (int N = base; N < cap; N += increment)
			{
				wave = waveGen.sinwave(a, numSamples, c, d);
				dftA = bdft.getDft(wave);
				dftB = dft4.getDft(wave);
				assertSimilar(dftA, dftB, error);

				wave = waveGen.coswave(a, numSamples, c, d);
				dftA = bdft.getDft(wave);
				dftB = dft4.getDft(wave);
				assertSimilar(dftA, dftB, error);

				a += 7.;
				c += 0.009;
				d += 0.05;
			}
		}
		TEST_METHOD(DFT5TEST)
		{
			DFT5 dft5 = DFT5();
			numSamples = 5;
			error = 0.000000001;
			int base = 1;
			int cap = 100;
			int increment = 1;
			double a = 10.;
			double c = 0.01;
			double d = 1.;
			for (int N = base; N < cap; N += increment)
			{
				wave = waveGen.sinwave(a, numSamples, c, d);
				dftA = bdft.getDft(wave);
				dftB = dft5.getDft(wave);
				assertSimilar(dftA, dftB, error);

				wave = waveGen.coswave(a, numSamples, c, d);
				dftA = bdft.getDft(wave);
				dftB = dft5.getDft(wave);
				assertSimilar(dftA, dftB, error);

				a += 7.;
				c += 0.009;
				d += 0.05;
			}
		}
		TEST_METHOD(DFT7TEST)
		{
			DFT7 dft7 = DFT7();
			numSamples = 7;
			error = 0.000000001;
			int base = 1;
			int cap = 100;
			int increment = 1;
			double a = 10.;
			double c = 0.01;
			double d = 1.;
			for (int N = base; N < cap; N += increment)
			{
				wave = waveGen.sinwave(a, numSamples, c, d);
				dftA = bdft.getDft(wave);
				dftB = dft7.getDft(wave);
				assertSimilar(dftA, dftB, error);

				wave = waveGen.coswave(a, numSamples, c, d);
				dftA = bdft.getDft(wave);
				dftB = dft7.getDft(wave);
				assertSimilar(dftA, dftB, error);

				a += 7.;
				c += 0.009;
				d += 0.05;
			}
		}
	};

	TEST_CLASS(PFFTTests)
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
		TEST_METHOD(SinWaveTest68)
		{
			error = 0.000000001;
			wave = waveGen.sinwave(10., 68, 0.01, 1.);
			dftA = bdft.getDft(wave);
			dftB = pfft.getDft(wave);
			assertSimilar(dftA, dftB, error);
		}
		TEST_METHOD(SinWaveTest82)
		{
			error = 0.000000001;
			wave = waveGen.sinwave(10., 82, 0.01, 1.);
			dftA = bdft.getDft(wave);
			dftB = pfft.getDft(wave);
			assertSimilar(dftA, dftB, error);
		}
		TEST_METHOD(SinWaveTest145)
		{
			error = 0.000000001;
			wave = waveGen.sinwave(10., 145, 0.01, 1.);
			dftA = bdft.getDft(wave);
			dftB = pfft.getDft(wave);
			assertSimilar(dftA, dftB, error);
		}
		TEST_METHOD(StressTest1)
		{
			error = 0.000000001;
			int base = 1;
			int cap = 100;
			int increment = 1;
			for (int N = base; N < cap; N += increment)
			{
				wave = waveGen.coswave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
			}
		}
		TEST_METHOD(StressTestXtreme)
		{
			error = 0.000000001;
			int base = 1;
			int cap = 100;
			int increment = 1;
			double a = 10.;
			double c = 0.01;
			double d = 1.;
			for (int N = base; N < cap; N += increment)
			{
				wave = waveGen.sinwave(a, N, c, d);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);

				wave = waveGen.coswave(a, N, c, d);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);

				a += 7.;
				c += 0.009;
				d += 0.05;
			}
		}
	};

	TEST_CLASS(CooleyTukeyTests)
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
		Logger logger = Logger();
		std::stringstream strStream = std::stringstream();
		std::string stringBuffer;
		WaveGenerator waveGen = WaveGenerator();
		CooleyTukeyFFT ctfft = CooleyTukeyFFT();
		BruteForceDft bdft = BruteForceDft();
		std::vector<std::complex<double>> dftA;
		std::vector<std::complex<double>> dftB;
		std::vector<std::complex<double>> wave;
		double error;
		TEST_METHOD(CosWaveTest32)
		{
			error = 0.000000001;
			wave = waveGen.coswave(10., 32, 0.01, 1.);
			dftA = bdft.getDft(wave);
			dftB = ctfft.getDft(wave);
			assertSimilar(dftA, dftB, error);
		}
		TEST_METHOD(SinWaveTest32)
		{
			error = 0.000000001;
			wave = waveGen.sinwave(10., 32, 0.01, 1.);
			dftA = bdft.getDft(wave);
			dftB = ctfft.getDft(wave);
			assertSimilar(dftA, dftB, error);
		}
		TEST_METHOD(CosWaveTest16)
		{
			error = 0.000000001;
			wave = waveGen.coswave(10., 16, 0.01, 1.);
			dftA = bdft.getDft(wave);
			dftB = ctfft.getDft(wave);
			assertSimilar(dftA, dftB, error);
		}
		TEST_METHOD(SinWaveTest16)
		{
			error = 0.000000001;
			wave = waveGen.sinwave(10., 16, 0.01, 1.);
			dftA = bdft.getDft(wave);
			dftB = ctfft.getDft(wave);
			assertSimilar(dftA, dftB, error);
		}
		TEST_METHOD(StressTest1)
		{
			error = 0.01;
			int base = 2;
			int cap = 12;
			int increment = 1;
			int N = 2;

			for (int i = base; i < cap; i += increment)
			{
				wave = waveGen.coswave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = ctfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
				N *= 2;
				strStream << "Status: " << i << " of " << cap << "...";
				stringBuffer = strStream.str();
				logger.WriteMessage(stringBuffer.c_str());
				strStream.str("");
				strStream.clear();
			}
		}
		TEST_METHOD(StressTestXtreme)
		{
			error = 0.000000001;
			int base = 2;
			int cap = 10;
			int increment = 1;

			int subBase = 1;
			int subCap = 100;
			int subIncrement = 1;

			int counter = 0;
			int numTests = cap * subCap;

			double a = 10.;
			double c = 0.01;
			double d = 1.;
			int N = 2;
			for (int i = base; i < cap; i += increment)
			{
				for (int ii = subBase; ii < subCap; ii += subIncrement)
				{
					wave = waveGen.sinwave(a, N, c, d);
					dftA = bdft.getDft(wave);
					dftB = ctfft.getDft(wave);
					assertSimilar(dftA, dftB, error);

					wave = waveGen.coswave(a, N, c, d);
					dftA = bdft.getDft(wave);
					dftB = ctfft.getDft(wave);
					assertSimilar(dftA, dftB, error);

					a += 7.;
					c += 0.009;
					d += 0.05;


					counter++;
				}
					strStream << "Status: " << i << " of " << cap << "...";
					stringBuffer = strStream.str();
					logger.WriteMessage(stringBuffer.c_str());
					strStream.str("");
					strStream.clear();
				N *= 2;
			}
		}
	};

	TEST_CLASS(PFFTStressTests)
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
		TEST_METHOD(StressTestPow3)
		{
			error = 0.000000001;
			int base = 3;
			int cap = 2000;
			int mult = 3;
			for (int N = base; N < cap; N *= mult)
			{
				wave = waveGen.sinwave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
			}
		}
		TEST_METHOD(StressTestPow2)
		{
			error = 0.000000001;
			int base = 2;
			int cap = 1000;
			int mult = 2;
			for (int N = base; N < cap; N *= mult)
			{
				wave = waveGen.sinwave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
			}
		}
		TEST_METHOD(StressTestPow5)
		{
			error = 0.000000001;
			int base = 5;
			int cap = 1000;
			int mult = 5;
			for (int N = base; N < cap; N *= mult)
			{
				wave = waveGen.sinwave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
			}
		}
		TEST_METHOD(StressTestPow7)
		{
			error = 0.000000001;
			int base = 7;
			int cap = 1000;
			int mult = 7;
			for (int N = base; N < cap; N *= mult)
			{
				wave = waveGen.sinwave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
			}
		}
		TEST_METHOD(StressTestPow6)
		{
			error = 0.000000001;
			int base = 6;
			int cap = 3000;
			int mult = 6;
			for (int N = base; N < cap; N *= mult)
			{
				wave = waveGen.sinwave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
			}
		}
		TEST_METHOD(StressTestMult2From3)
		{
			error = 0.000000001;
			int base = 3;
			int cap = 3000;
			int mult = 2;
			for (int N = base; N < cap; N *= mult)
			{
				wave = waveGen.sinwave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
			}
		}
		TEST_METHOD(StressTestMult3From2)
		{
			error = 0.000000001;
			int base = 2;
			int cap = 3000;
			int mult = 3;
			for (int N = base; N < cap; N *= mult)
			{
				wave = waveGen.sinwave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
			}
		}
		TEST_METHOD(StressTestMult2From5)
		{
			error = 0.000000001;
			int base = 5;
			int cap = 3000;
			int mult = 2;
			for (int N = base; N < cap; N *= mult)
			{
				wave = waveGen.sinwave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
			}
		}
		TEST_METHOD(StressTestMult2From7)
		{
			error = 0.000000001;
			int base = 7;
			int cap = 3000;
			int mult = 2;
			for (int N = base; N < cap; N *= mult)
			{
				wave = waveGen.sinwave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
			}
		}
		TEST_METHOD(StressTestMult5From2)
		{
			error = 0.000000001;
			int base = 2;
			int cap = 3000;
			int mult = 5;
			for (int N = base; N < cap; N *= mult)
			{
				wave = waveGen.sinwave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
			}
		}
		TEST_METHOD(StressTestMult7From2)
		{
			error = 0.000000001;
			int base = 2;
			int cap = 3000;
			int mult = 7;
			for (int N = base; N < cap; N *= mult)
			{
				wave = waveGen.sinwave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
			}
		}
		TEST_METHOD(StressTestMult5From3)
		{
			error = 0.000000001;
			int base = 3;
			int cap = 3000;
			int mult = 5;
			for (int N = base; N < cap; N *= mult)
			{
				wave = waveGen.sinwave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
			}
		}
		TEST_METHOD(StressTestMult3From5)
		{
			error = 0.000000001;
			int base = 5;
			int cap = 3000;
			int mult = 3;
			for (int N = base; N < cap; N *= mult)
			{
				wave = waveGen.sinwave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
			}
		}
		TEST_METHOD(StressTestMult7From3)
		{
			error = 0.000000001;
			int base = 3;
			int cap = 3000;
			int mult = 7;
			for (int N = base; N < cap; N *= mult)
			{
				wave = waveGen.sinwave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
			}
		}
		TEST_METHOD(StressTestMult3From7)
		{
			error = 0.000000001;
			int base = 7;
			int cap = 3000;
			int mult = 3;
			for (int N = base; N < cap; N *= mult)
			{
				wave = waveGen.sinwave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
				wave = waveGen.coswave(10., N, 0.01, 1.);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
			}
		}
		TEST_METHOD(StressTestMult3From7Xtreme)
		{
			error = 0.000000001;
			int base = 7;
			int cap = 3000;
			int mult = 3;
			double a = 10.;
			double c = 0.01;
			double d = 1.;
			for (int N = base; N < cap; N *= mult)
			{
				wave = waveGen.sinwave(a, N, c, d);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
				wave = waveGen.coswave(a, N, c, d);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);

				a += 7.;
				c += 0.009;
				d += 0.05;
			}
		}
		TEST_METHOD(StressTestLarge)
		{
			error = 0.000000001;
			int base = 100;
			int cap = 300;
			int increment = 8;
			double a = 10.;
			double c = 0.01;
			double d = 1.;
			for (int N = base; N < cap; N += increment)
			{
				wave = waveGen.sinwave(a, N, c, d);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);
				wave = waveGen.coswave(a, N, c, d);
				dftA = bdft.getDft(wave);
				dftB = pfft.getDft(wave);
				assertSimilar(dftA, dftB, error);

				a += 1.;
				c += 0.001;
				d += 0.01;
			}
		}
	};
}