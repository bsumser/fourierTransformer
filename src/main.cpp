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
#include <unistd.h>
#include <fstream>
#include <pthread.h>
#include <omp.h>
#include "../include/timelord.h"
#include "../include/wavio.h"
#include "../include/imageLoader.h"


using namespace std::chrono;
using namespace std;


vector<float> floatGenerator(int vectorSize);
int processArgs(int argc, char* argv[]);
void writeToFile(string file, vector<complex<double>> output);
void writeVectorToFile(string file, vector<double> output);
vector<complex<double>> discreteFourierTransform(vector<double> input);
vector<complex<double>> parallelDiscreteFourierTransform(vector<double> input);
vector<complex<double>> discreteFourierTransformTurkey(vector<double> input);
vector<double> signalGenerator(int sampleSize);
vector<double> detectPitches(vector<complex<double>> output);
unsigned int bitReverse(unsigned int input, int n);


int n = 1000;
double test = 0.000000001;
bool wav = false;
bool image = false;
bool num = false;
const char* inFile = "wavs/a.wav";


int main(int argc, char* argv[]) {
    if (processArgs(argc, argv))
		return 1;

    //get path for image
    char* path = argv[2];

    //attempting to load image
    ImageLoader imageLoader(path);
    imageLoader.grayscaler();
    discreteCosineTransform(imageLoader);

    // Init Variables
    high_resolution_clock::time_point start, stop;
    Timelord newTimeLord();
    duration<double> duration;
    vector<double> input;
	vector<complex<double>> output;

	if (wav) {
		cout << "Reading " << inFile << ":" << endl;
		start = high_resolution_clock::now();
		input = read_wav(inFile);
		stop = high_resolution_clock::now();
		cout << "done in " << duration.count() << " seconds" << endl;
		n = input.size();
		writeVectorToFile("sig.txt", input);
	}

	if (image) {
		cout << "Looking at " << inFile << " ... ";
		start = high_resolution_clock::now();
		stop = high_resolution_clock::now();
		cout << "done in " << duration.count() << " seconds" << endl;
		n = input.size();
	}

	cout << "N = " << n << endl;	
	
	if (!wav && !image) {
	    // Signal Generator
		cout << "Generating signal ... ";
	    start = high_resolution_clock::now();
	    input = signalGenerator(n);
	    stop = high_resolution_clock::now();
	    duration = stop - start;
	    cout << "done in " << duration.count() << " seconds" << endl << endl;
		writeVectorToFile("sig.txt", input);
	}

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


    //DFTCT
    cout << "Starting DFT Cooley-Turkey ... ";
    start = high_resolution_clock::now();
    output = discreteFourierTransformTurkey(input);
    stop = high_resolution_clock::now();
    duration = stop - start;
    cout << "done in " << duration.count() << " seconds" << endl;
    writeToFile("ct.txt", output);
   
	if (wav) {
		// Pitch Detection
		cout << "Detecting Pitches ... ";
		start = high_resolution_clock::now();
		vector<double> pitches = detectPitches(output);
		stop = high_resolution_clock::now();
		duration = stop - start;
		cout << "done in " << duration.count() << " seconds" << endl;
		cout << "Pitches detected:" << endl;
		for (unsigned int i = 0; i < pitches.size(); i++)
			cout << pitches[i] << endl;
	}


	return 0;
}

int processArgs(int argc, char* argv[])
{
	int opt;
	while ((opt = getopt(argc, argv, "w::i::n:")) != -1) {
		switch (opt) {
			case 'w':
				wav = true;
				if (optarg == NULL && optind < argc && argv[optind][0] != '-')
					optarg = argv[optind++];
				if (optarg != NULL)
					inFile = optarg;
				break;
			case 'i':
				image = true;
				if (optarg == NULL && optind < argc && argv[optind][0] != '-')
					optarg = argv[optind++];
				if (optarg != NULL)
					inFile = optarg;
				break;
			case 'n':
				num = true;	
				n = stoi(optarg);
				break;
			case '?':
				return 1;
		}
	}
	return 0;
}

void writeToFile(string file, vector<complex<double>> vec)
{
    #ifdef PRINT
		cout << "Writing to file ... ";
		ofstream f;
		f.open("out/" + file, ios::trunc);
		for (int i = 0 ; i < n; i++)
			f << vec[i].real() << ", " << vec[i].imag() << endl;
		f.close();
		cout << "done" << endl << endl;
    #endif
}

void writeVectorToFile(string file, vector<double> vec)
{
    #ifdef PRINT
		cout << "Writing to file ... ";
		ofstream f;
		f.open("out/" + file, ios::trunc);
		for (int i = 0 ; i < n; i++)
			f << vec[i] << endl;
		f.close();
		cout << "done" << endl << endl;
    #endif
}

void discreteCosineTransform(ImageLoader imageloader)
{
    int length = imageloader.image.size();
    int coeff = (2/length)^(0.5);
    for (int i = 0; i < length; i++) {
        //For each pixel in image, apply the discrete cosine transform
        //to access pixels use (int)imageloader.image[i]
        //the image data is a 1d vector of char style pixels

    }
    //for (int i = 0; i < imageloader.image.size(); i+=4)
    //{
    //    std::cout << "Pixel " << i << " is : " 
    //    << "r:" << (int)imageloader.image[i] << " " 
    //    << "g:" << (int)imageloader.image[i + 1] << " "
    //    << "b:" << (int)imageloader.image[i + 2] << " "
    //    << "a:" << (int)imageloader.image[i + 3] << " "
    //    << std::endl;
    //}
}

vector<complex<double>> discreteFourierTransform(vector<double> input)
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

vector<complex<double>> parallelDiscreteFourierTransform(vector<double> input)
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

vector<complex<double>> discreteFourierTransformTurkey(vector<double> input)
{
    vector<complex<double>> output;

	int high = ceil(log(input.size()) / log(2));
	int pad = pow(2, high);

	while ((int)input.size() < pad)
		input.push_back(0.0);
    
	unsigned int N = input.size();
    output.reserve(N);

	for (unsigned int i = 0; i < N; i++) {
		int reversed = bitReverse(i, pad);
		output.push_back(input[reversed]);
	}

	cout << N << endl;

	for (int s = 1; s <= pad; s++) {
		cout << "s: " << s << endl;
		int m = pow(2, s);
		double real = cos((-2 * M_PI) / m);
		double imag = -sin((-2 * M_PI) / m);
		complex<double> wm;
		wm.real(real);
		wm.imag(imag);
		for (int k = 0; k < pad; k++) {
			cout << "k: " << k << endl;
			complex<double> w;
			w.real(1.0);
			w.imag(0.0);
			for (int j = 0; j < (m / 2) - 1; j++) {
				cout << "j: " << j << endl;
				complex<double> t = w * output[k + j + m / 2];
				complex<double> u = output[k + j];
				output[k + j] = u + t;
				output[k + j + m / 2] = u - t;
				w *= wm;
			}
		}
	}
    return output;
}

unsigned int bitReverse(unsigned int num, int N)
{
	int n = 0;
	for (int i = 0; i < N; i++) {
		n <<= 1;
		n |= (num & 1);
		num >>= 1;
	}
	return n;
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

vector<double> signalGenerator(int sampleSize)
{
    double step = 1.0 / (double)sampleSize;
	double test = 0.00000000001;

    vector<double> output; //output vector to store the test signal
    output.reserve(sampleSize);  //allocate proper size for vector

    //#pragma omp parallel for
    for (double i = 0.0; i < (double)1.0; i += step) {
		double tmp = sin((double)880.0 * (double)M_PI * i);
		if (fabs(tmp) < test) 
			tmp = 0.0;
        output.push_back(tmp);
    }
    
    return output;
}

vector<double> detectPitches(vector<complex<double>> output)
{
	vector<double> pitches;
	for (int i = 0; i < (int)output.size(); i++) {
		if (fabs(output[i].imag() > 500.0))
			pitches.push_back(i);
	}
	return pitches;
}
