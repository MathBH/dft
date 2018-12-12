

	Tesla Belzile-Ha														December 5th 2018
	----------------														-----------------


							   :Fast Fourier Transform Algorithms:
							   -\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-


		Implementations of Fast Fourier Transform Algorithms:
		=====================================================

				Currently includes:
				-------------------
							0) Brute-Force DFT Algorithm
							1) Cooley-Tukey FFT Algorithm
							2) Prime-Factor FFT Algorithm

--------------------------------------------------------------------------------------------------
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

 
[0] Brute-Force DFT Algorithm:
==============================

	A class built from the "DFTAlgorithm" template. Calculates a sample set's Discrete Fourier
	Transform using a brute-force implementation.

	Brute Force Algorithm based on equations from "Numerical Recipes in C++" p.508.
				

				Files:
				======

						"BruteForceDft.h", "BruteForceDft.cpp"


				Functions:
				==========
						
						reverseDft(dft):
						----------------

							Rebuilts the original signal from the input dft.

							input:
									vector<complex<double>> dft:

												discrete fourier transform of a signal

							output:
									vector<complex<double>> signal:

												signal constructed from the input dft
				Templates: 
				=========

						"DFTAlgorithm.h":
						-----------------

						Template denoting a class that can calculate the discrete fourier
						transform for a given set of samples


						Functions:
						==========
						
							  getDft(samples):
							  ----------------

								Calculates the discrete fourier transform of given samples.


								input:
										  vector<complex<double>> samples:
										
													 Data set from which the dft is to be
													 calculated.
										
								output:

										 vector<complex<double>> dft:

													Discrete fourier transform calculated
													from the input samples.



 
[1] Cooley-Tukey FFT Algorithm:
===============================

	Class with the "DFTAlgorithm" template. Code modified from "Numerical Recipes in C++" 
	p.513. Calculates a sample set's Discrete Fourier Transform using the Cooley-Tukey Algorithm.

				Files:
				======

						"CooleyTukeyFFT.h", "CooleyTukeyFFT.cpp"


				Templates: 
				=========

						"DFTAlgorithm.h":
						-----------------

						Template denoting a class that can calculate the discrete fourier
						transform for a given set of samples


						Functions:
						==========

							  getDft(samples):
							  ----------------

								Calculates the discrete fourier transform of given samples.


								input:
										  vector<complex<double>> samples:
										
													 Data set from which the dft is to be
													 calculated.
										
								output:

										 vector<complex<double>> dft:

													Discrete fourier transform calculated
													from the input samples.


 

[2] Prime-Factor FFT Algorithm:
===============================

	Class with the "DFTAlgorithm" template. Calculates a sample set's Discrete Fourier
	Transform using the Prime-Factor algorithm.

	Algorithm from "Fast Fourier Transform and Convolution Algorithms" p.125-129.

				Files:
				======

						"PrimeFactorDft.h", "PrimeFactorDft.cpp"


				Templates:
				=========

						"DFTAlgorithm.h":
						-----------------

						Template denoting a class that can calculate the discrete fourier
						transform for a given set of samples


						Functions:
						==========
						
							  getDft(samples):
							  ----------------

								Calculates the discrete fourier transform of given samples.


								input:
										  vector<complex<double>> samples:
										
													 Data set from which the dft is to be
													 calculated.
										
								output:

										 vector<complex<double>> dft:

													Discrete fourier transform calculated
													from the input samples.


				Related Classes:
				================
					
						The Prime-Factor Algorithm splits up the data into smaller data sets
						which it then passes to sub-algorithms. The sub algorithms used are
						as follows,

						- BruteForceDft:

								used as last resort

						- CooleyTukeyFFT:

								used when data subset is of a size that does not match any other
								subalgorithm but is a power of 2.

						- <SmallDFTs>:

								Set of classes for calculating small dfts. Algorithms from
								"Fast Fourier Transform and Convolution" p.145 - p150.
								These are used by the Prime-Factor Algorithm along with
								the Brute-Force Algorithm and Cooley-Tukey Algorithm.
								
								When splitting the data into factors, any factor that matches
								the small dfts is calculated with them, any factor that does
								not match but is a power of two is calculated using the
								Cooley-Tukey Algorithm. All other factor's corresponding
								subsets are calculated using the Brute-Force algorithm.

								Files:
										"SmallDFTs.h", "SmallDFTs.cpp"

								Classes:
									- DFT2: 2-point dft evaluator - template: DFTAlgorithm
									- DFT3: 3-point dft evaluator - template: DFTAlgorithm
									- DFT4:	4-point dft evaluator - template: DFTAlgorithm
									- DFT5:	5-point dft evaluator - template: DFTAlgorithm
									- DFT7:	7-point dft evaluator - template: DFTAlgorithm


--------------------------------------------------------------------------------------------------
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\