/* Program driver for a series of Parallel Fourier Tranformations.
 *
 * Authors:
 *	Brett Sumser
 *	Justin Spidell
 */
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <chrono>
#include <stdio.h>
#include <string>
#include <fstream>
#include <pthread.h>
#include <omp.h>
#include "../include/timelord.h"
#include "../include/imageLoader.h"


using namespace std::chrono;
using namespace std;


vector<float> floatGenerator(int vectorSize);
void processArgs(int argc, char* argv[]);
void writeToFile(string file, vector<complex<double>> output);
void discreteCosineTransform(ImageLoader imageloader);
vector<complex<double>> discreteFourierTransform(vector<complex<double>> input);
vector<complex<double>> parallelDiscreteFourierTransform(vector<complex<double>> input);
vector<complex<double>> discreteFourierTransformTurkey(vector<complex<double>> input);
vector<complex<double>> signalGenerator(int sampleSize);


int n = 1000;
double test = 0.000000001;


int main(int argc, char* argv[]) {
    // Get N
    processArgs(argc, argv);
    cout << "N = " << n << endl;

    //get path for image
    char* path = argv[2];

    //attempting to load image
    ImageLoader imageLoader(path);
    imageLoader.grayscaler();
    discreteCosineTransform(imageLoader);

    // Init Variables
    high_resolution_clock::time_point start, stop;
    duration<double> duration;
    vector<complex<double>> input, output;
    Timelord newTimeLord();

    // Signal Generator
    cout << "Generating signal ... ";
    start = high_resolution_clock::now();
    input = signalGenerator(n);
    stop = high_resolution_clock::now();
    duration = stop - start;
    cout << "done in " << duration.count() << " seconds" << endl << endl;
	writeToFile("sig.txt", input);


    // DFT
    cout << "Starting DFT ... ";
    start = high_resolution_clock::now();
    output = discreteFourierTransform(input);
    stop = high_resolution_clock::now();
    duration = stop - start;
    cout << "done in " << duration.count() << " seconds" << endl;
    writeToFile("dft.txt", output);

    // DFTF
    cout << "Starting // DFT ... ";
    start = high_resolution_clock::now();
    output = parallelDiscreteFourierTransform(input);
    stop = high_resolution_clock::now();
    duration = stop - start;
    cout << "done in " << duration.count() << " seconds" << endl;
    writeToFile("pdft.txt", output);


    //DFTT
    cout << "Starting DFT Cooley-Turkey ... ";
    start = high_resolution_clock::now();
    output = discreteFourierTransformTurkey(input);
    stop = high_resolution_clock::now();
    duration = stop - start;
    cout << "done in " << duration.count() << " seconds" << endl;
    writeToFile("ct.txt", output);


    return 0;
}

void processArgs(int argc, char* argv[])
{
    // TODO: Add a flag reader
    if (argc > 1) {
	//n = stoi(argv[1]);
    for (int i = 0; i < argc; i++) {
        std::cout << argv[i] << std::endl;
    }
    }
}

void writeToFile(string file, vector<complex<double>> output)
{
    #ifdef PRINT
		cout << "Writing to file ... ";
		ofstream f;
		f.open("out/" + file, ios::trunc);
		for (int i = 0 ; i < n; i++)
			f << output[i].real() << ", " << output[i].imag() << endl;
		f.close();
		cout << "done" << endl << endl;
    #endif
}

void discreteCosineTransform(ImageLoader imageloader)
{
    for (int i = 0; i < imageloader.image.size(); i+=4)
    {
        std::cout << "Pixel " << i << " is : " 
        << "r:" << (int)imageloader.image[i] << " " 
        << "g:" << (int)imageloader.image[i + 1] << " "
        << "b:" << (int)imageloader.image[i + 2] << " "
        << "a:" << (int)imageloader.image[i + 3] << " "
        << std::endl;
    }
}

vector<complex<double>> discreteFourierTransform(vector<complex<double>> input)
{
    //initialize sizes of samples
    int N = input.size();
    int K = input.size();

	int k;
	int n;


    //init variable for internal loop
    complex<double> innerSum = complex<double>(0.0,0.0);

    //init vector of std::complex doubles for results
    vector<complex<double>> output;

    //set output vector to have enough space
    output.reserve(N);

	complex<double> tmp = complex<double>(0.0,0.0);

	double real = 0.0;
	double imag = 0.0;

	//sigma notation algorithm start
    for (k = 0; k < K; k++) {
        innerSum = complex<double>(0.0,0.0);
		for (n = 0; n < N; n++) {
	        real = cos((2 * M_PI * k * n) / N);
            imag = -sin((2 * M_PI * k * n) / N);

			tmp.real(real);
			tmp.imag(imag);
			tmp *= input[n];
			innerSum += tmp;
		}
		if (fabs(innerSum.real()) < test)
			innerSum.real(0.0);
		if (fabs(innerSum.imag()) < test)
			innerSum.imag(0.0);
	    output.push_back(innerSum);
    }
    
    return output;
}

vector<complex<double>> parallelDiscreteFourierTransform(vector<complex<double>> input)
{
    //initiliaze sizes of samples
    int N = input.size();
    int K = input.size();

    //loop indexes
    int k;
    int n;

    //init variable for internal loop
    complex<double> innerSum = complex<double>(0.0,0.0);

    //init vector of std::complex doubles for results
    vector<complex<double>> output;

    //set output vector to have enough space
    output.resize(N);

	complex<double> tmp = complex<double>(0.0,0.0);

	double real = 0.0;
	double imag = 0.0;

	#pragma omp parallel for default(none) private(real, imag, innerSum, k, tmp, n) shared(input, output, K, N, cout, test)
	for (k = 0; k < K; k++) {
		innerSum = complex<double>(0.0,0.0);	
		for (n = 0; n < N; n++) {
			real = cos(((2 * M_PI) / N) * k * n);
			imag = -sin(((2 * M_PI) / N) * k * n);

			tmp.real(real);
			tmp.imag(imag);
			tmp *= input[n];
			innerSum += tmp;
		}
		if (fabs(innerSum.real()) < test)
			innerSum.real(0.0);
		if (fabs(innerSum.imag()) < test)
			innerSum.imag(0.0);
		#pragma omp critical
		output.at(k) = innerSum;
	}

    return output;
}

vector<complex<double>> discreteFourierTransformTurkey(vector<complex<double>> input)
{
    //initiliaze sizes of samples
    int N = input.size();
    int K = input.size();

    //init variable for internal loop
    complex<double> innerSum;

    //init vector of std::complex doubles for results
    vector<complex<double>> output;

    //set output vector to have enough space
    output.reserve(N);
    
    for (int k = 0; k < K; k++)
    {
        innerSum = complex<double>(0,0);
        for (int n = 0; n < N / 2 - 1; n++) {
            // process real and imaginary parts of sum via definition in
            // "https://en.wikipedia.org/wiki/Discrete_Fourier_transform"
            double real = cos(((2 * M_PI) / (N / 2)) * k * (2 * n));
            double imag = sin(((2 * M_PI) / (N / 2)) * k * (2 * n));

            //store as std::complex double for multiplication
            complex<double> w (real, -imag);
            innerSum += input[n] * w;
        }
        for (int n = 1; n < N / 2 - 1; n++) {
            // process real and imaginary parts of sum via definition in
            // "https://en.wikipedia.org/wiki/Discrete_Fourier_transform"
            double real = cos(((2 * M_PI) / (N / 2)) * k * (2 * n + 1));
            double imag = sin(((2 * M_PI) / (N / 2)) * k * (2 * n + 1));

            //store as std::complex double for multiplication
            complex<double> w (real, -imag);
            innerSum += input[n] * w;
        }
        //add value to the output vector
        output.push_back(innerSum);
    }
    
    return output;
}

vector<float> floatGenerator(int vectorSize)
{
    //initialize stuff fo random float generator
    default_random_engine floatGenerator;
    mt19937 mt(floatGenerator());
    uniform_real_distribution<float> uniform_distance(0, 10.001);
    
    //float vector to return
    vector<float> floatVector;

    //#pragma omp parallel for 
    for (int i = 0; i < vectorSize - 1; i++) {
        floatVector.push_back(uniform_distance(mt));
        //cout << "rand is " << uniform_distance(mt) << endl;
    }
    return floatVector;
}

vector<complex<double>> signalGenerator(int sampleSize)
{
    int N = sampleSize;
    double amp = 3;  //amplitude, peak deviation from 0
    double freq = 0;   //ordinary frequency, oscillations
    double time = 0;   //time
    double shift = M_PI / 2;  //phase shift of signal

	double test = 0.00000000001;

    vector<complex<double>> output; //output vector to store the test signal
    output.reserve(N);  //allocate proper size for vector

    //#pragma omp parallel for
    for (int i = 0; i < N; i++) {
		double tmp = amp * cos((2 * M_PI / N) * i + shift);
		if (fabs(tmp) < test) 
			tmp = 0.0;
		auto sample = complex<double>(tmp, 0.0);
        output.push_back(sample);
    }
    
    return output;
}
