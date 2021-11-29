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


int processArgs(int argc, char* argv[]);
void writeToFile(string file, vector<complex<double>> output);
void writeVectorToFile(string file, vector<double> output);
vector<complex<double>> discreteFourierTransform(vector<double> input);
vector<complex<double>> parallelDiscreteFourierTransform(vector<double> input);
vector<complex<double>> discreteFourierTransformTurkey(vector<double> input);
vector<double> signalGenerator(int sampleSize);
vector<double> detectPitches(vector<complex<double>> output);
unsigned int bitReverse(unsigned int input, int log2n);
void discreteCosineTransform(ImageLoader imageloader);


int n = 1024;
double test = 0.000000001;
bool wav = false;
bool image = false;
bool num = false;
bool err = false;
string inFile = "wavs/a.wav";


int main(int argc, char* argv[]) {
    if (processArgs(argc, argv))
		return 1;

    // Init Variables
    high_resolution_clock::time_point start, stop;
    Timelord newTimeLord();
    duration<double> duration;
    vector<double> input;
	vector<complex<double>> dftout, pdftout, dftctout, dctout;

	if (wav) {
		cout << "Reading " << inFile << ":" << endl;
		start = high_resolution_clock::now();
		input = read_wav(inFile.c_str());
		stop = high_resolution_clock::now();
		cout << "done in " << duration.count() << " seconds" << endl;
		n = input.size();
		writeVectorToFile("sig.txt", input);
	}

	if (image) {
		cout << "Looking at " << inFile << " ... ";
		//attempting to load image
		ImageLoader imageLoader(inFile.c_str());
		imageLoader.grayscaler();
		start = high_resolution_clock::now();
		discreteCosineTransform(imageLoader);
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
    dftout = discreteFourierTransform(input);
    stop = high_resolution_clock::now();
    duration = stop - start;
    cout << "done in " << duration.count() << " seconds" << endl;
    writeToFile("dft.txt", dftout);


    // PDFT
    cout << "Starting // DFT ... ";
    start = high_resolution_clock::now();
    pdftout = parallelDiscreteFourierTransform(input);
    stop = high_resolution_clock::now();
    duration = stop - start;
    cout << "done in " << duration.count() << " seconds" << endl;
    writeToFile("pdft.txt", pdftout);


    //DFTCT
    cout << "Starting DFT Cooley-Turkey ... ";
    start = high_resolution_clock::now();
    dftctout = discreteFourierTransformTurkey(input);
    stop = high_resolution_clock::now();
    duration = stop - start;
    cout << "done in " << duration.count() << " seconds" << endl;
    if (!err)
		writeToFile("ct.txt", dftctout);
	else
		cout << endl;

	if (wav) {
		// Pitch Detection
		cout << "Detecting Pitches ... ";
		start = high_resolution_clock::now();
		vector<double> pitches = detectPitches(dftout);
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
	if (wav && image) {
		cout << "I can't do both an image and a wav file... sorry" << endl;
		return 1;
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
    // pow? org: (2/length)^(0.5)
	//int coeff = pow((2 / length), (0.5));
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

	int high = log(input.size()) / log(2);
    
	double tmp = log(input.size()) / log(2);

	if (high != tmp) {
		cout << endl << "This algorithm only works with power of 2 sized arrays!" << endl;
		err = true;
		return output;
	}

	unsigned int N = (unsigned int)input.size();
    output.reserve(N);

	int top = input.size();

	for (unsigned int i = 0; i < N; ++i) {
		unsigned int reversed = bitReverse(i, high);
		complex<double> tmp;
		tmp.real(input[reversed]);
		tmp.imag(0.0);
		output.push_back(tmp);
	}

	const complex<double> iota (0, 1);
	for (int s = 1; s <= high; ++s) {
		int m = 1 << s;
		int m2 = m >> 1;
		complex<double> w (1.0, 0.0);
		double real = cos(M_PI / m2);
		double imag = -sin(M_PI / m2);
		complex<double> wm (real, imag);
		if (fabs(wm.real()) < test)
			wm.real(0.0);
		if (fabs(wm.imag()) < test)
			wm.imag(0.0);
		for (int j = 0; j < m2; ++j) {
			complex<double> t, u;
			int k;
			// #pragma omp parallel for default(none) shared(j, output, test, top, m, m2, w) private(t, u, k)	
			for (k = j; k < top; k += m) {
				t = (w * output[k + m2]);
				u = output[k];
				// #pragma omp critical
				{
				output[k] = u + t;
				output[k + m2] = u - t;
				if (fabs(output[k].real()) < test)
					output[k].real(0.0);
				if (fabs(output[k].imag()) < test)
					output[k].imag(0.0);
				if (fabs(output[k+m2].real()) < test)
					output[k+m2].real(0.0);
				if (fabs(output[k+m2].imag()) < test)
					output[k+m2].imag(0.0);
				}
			}
			w *= wm;
		}
	}
    return output;
}

unsigned int bitReverse(unsigned int num, int log2n)
{
	unsigned int rev = 0;

	for (int i = 0; i < log2n; i++) {
		rev <<= 1;
		rev |= (num & 1);
		num >>= 1;
	}

	return rev;
}

vector<double> signalGenerator(int sampleSize)
{
	double test = 0.00000000001;

    vector<double> output; //output vector to store the test signal
    output.resize(sampleSize);  //allocate proper size for vector

	double tmp;
	int i;

    #pragma omp parallel for default(none) shared(test, output, sampleSize) private(i, tmp)
    for (i = 0; i < sampleSize; i++) {
		tmp = sin(880.0 * M_PI * (double)((double)i / (double)sampleSize));
		if (fabs(tmp) < test) 
			tmp = 0.0;
        output.at(i) = tmp;
    }
    
    return output;
}

vector<double> detectPitches(vector<complex<double>> output)
{
	vector<double> pitches;
	for (int i = 0; i < (int)output.size() / 2; i++) {
		if (fabs(output[i].imag()) > 700.0)
			pitches.push_back(i);
	}
	return pitches;
}

