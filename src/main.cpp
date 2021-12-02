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
void writeVectorToFile(string file, vector<double> output, int N);
vector<complex<double>> DFT(vector<double> input);
vector<complex<double>> PDFT(vector<double> input);
vector<complex<double>> CT(vector<double> input);
vector<complex<double>> CTP1(vector<double> input);
vector<complex<double>> CTP2(vector<double> input);
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
vector<omp_lock_t> locks;


int main(int argc, char* argv[]) {
    if (processArgs(argc, argv))
		return 1;

    // Init Variables
    high_resolution_clock::time_point start, stop;
    Timelord newTimeLord();
    duration<double> duration;
	vector<double> durations;
    vector<double> input;
	vector<complex<double>> dftout, pdftout, ctout, ctp1out, ctp2out;

	if (wav) {
		cout << "Reading " << inFile << ":" << endl;
		start = high_resolution_clock::now();
		input = read_wav(inFile.c_str());
		stop = high_resolution_clock::now();
		cout << "done in " << duration.count() << " seconds" << endl;
		durations.push_back(duration.count());
		n = input.size();
		writeVectorToFile("sig.txt", input, n);
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
		durations.push_back(duration.count());
		n = input.size();
	}

	cout << "N = " << n << endl << endl;

	if (!wav && !image) {
	    // Signal Generator
		cout << "Generating signal ... ";
	    start = high_resolution_clock::now();
	    input = signalGenerator(n);
	    stop = high_resolution_clock::now();
	    duration = stop - start;
	    cout << "done in " << duration.count() << " seconds" << endl;
		durations.push_back(duration.count());
		writeVectorToFile("sig.txt", input, n);
	}

	// Locks
	cout << "Setting up locks ... ";
    start = high_resolution_clock::now();
	locks.reserve(n);
	for (int i = 0; i < n; i++)
		omp_init_lock(&(locks[i]));
    stop = high_resolution_clock::now();
	duration = stop - start;
	cout << "done in " << duration.count() << " seconds" << endl << endl;;
	durations.push_back(duration.count());

    // DFT
    cout << "Starting DFT ... ";
    start = high_resolution_clock::now();
    dftout = DFT(input);
    stop = high_resolution_clock::now();
    duration = stop - start;
    cout << "done in " << duration.count() << " seconds" << endl;
	durations.push_back(duration.count());
    writeToFile("dft.txt", dftout);


    // DFTP1
    cout << "Starting // DFT ... ";
    start = high_resolution_clock::now();
    pdftout = PDFT(input);
    stop = high_resolution_clock::now();
    duration = stop - start;
    cout << "done in " << duration.count() << " seconds" << endl;
	durations.push_back(duration.count());
    writeToFile("pdft.txt", pdftout);


    //CT
    cout << "Starting Cooley-Turkey ... ";
    start = high_resolution_clock::now();
    ctout = CT(input);
    stop = high_resolution_clock::now();
    duration = stop - start;
    cout << "done in " << duration.count() << " seconds" << endl;
	durations.push_back(duration.count());
	writeToFile("ct.txt", ctout);


    //CTP1
    cout << "Starting V1 // Cooley-Turkey ... ";
    start = high_resolution_clock::now();
    ctp1out = CTP1(input);
    stop = high_resolution_clock::now();
    duration = stop - start;
    cout << "done in " << duration.count() << " seconds" << endl;
	durations.push_back(duration.count());
	writeToFile("ctp1.txt", ctp1out);

	
    //CTP2
	cout << "Starting V2 // Cooley-Turkey ... ";
    start = high_resolution_clock::now();
    ctp2out = CTP2(input);
    stop = high_resolution_clock::now();
    duration = stop - start;
    cout << "done in " << duration.count() << " seconds" << endl;
	durations.push_back(duration.count());
	if (!err)
		writeToFile("ctp2.txt", ctp2out);

	for (int i = 0; i < n; i++)
		omp_destroy_lock(&(locks[i]));


	// MPI DFT

	// TODO MAKE MPI WORK WITH BASIC DFT
	// THEN COOLEY TUKEY, WILL HAVE TO PAD WITH 0'S IN MAIN
	// RATHER THAN IN THE FUNCTION
    MPI_Init(NULL, NULL);
	

    // Current rank's ID
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // Total number of ranks
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);


    // Read the sparse matrix and store it in row_ind, col_ind, and val,
    // also known as co-ordinate format (COO).
	vector<double> mpiInput = memcpy()

    // Rank 0 now determines how work will be distributed among the ranks
    int nnz_per_rank = 0;
    if(world_rank == 0) {
        nnz_per_rank = ( + world_size - 1) / world_size;
    }
    start = MPI_Wtime();
    // Broadcast this to everyone
    MPI_Bcast(&nnz_per_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Also broadcast m and n
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    end = MPI_Wtime();
    timer[VEC_BCAST_TIME] = end - start;


    // Now, let's send the sparse matrix
    // First, pad the data so that we can use MPI_Scatter instead of 
    // MPI_Scatterv
    if(world_rank == 0) {
        int new_nnz = nnz_per_rank * world_size;
        int* row_ind_tmp = (int*) malloc(sizeof(int) * new_nnz);
        assert(row_ind_tmp);
        memset(row_ind_tmp, 0, sizeof(int) * new_nnz);
        int* col_ind_tmp = (int*) malloc(sizeof(int) * new_nnz);
        assert(col_ind_tmp);
        memset(col_ind_tmp, 0, sizeof(int) * new_nnz);
        double* val_tmp = (double*) malloc(sizeof(double) * new_nnz);
        assert(val_tmp);
        memset(val_tmp, 0, sizeof(double) * new_nnz);

        memcpy(row_ind_tmp, row_ind, sizeof(int) * nnz);
        memcpy(col_ind_tmp, col_ind, sizeof(int) * nnz);
        memcpy(val_tmp, val, sizeof(double) * nnz);

        free(row_ind);
        free(col_ind);
        free(val);
        row_ind = row_ind_tmp;
        col_ind = col_ind_tmp;
        val = val_tmp;
    } else {
        // Everyone else should get ready to receive the appropriate 
        // amount of data
        row_ind = (int*) malloc(sizeof(int) * nnz_per_rank);
        assert(row_ind);
        col_ind = (int*) malloc(sizeof(int) * nnz_per_rank);
        assert(col_ind);
        val = (double*) malloc(sizeof(double) * nnz_per_rank);
        assert(val);
    }

    start = MPI_Wtime();    
    // Scatter the data to each node
    MPI_Scatter(row_ind, nnz_per_rank, MPI_INT, row_ind, nnz_per_rank, MPI_INT,
                0, MPI_COMM_WORLD);
    MPI_Scatter(col_ind, nnz_per_rank, MPI_INT, col_ind, nnz_per_rank, MPI_INT,
                0, MPI_COMM_WORLD);
    MPI_Scatter(val, nnz_per_rank, MPI_DOUBLE, val, nnz_per_rank, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    end = MPI_Wtime();
    timer[MAT_SCATTER_TIME] = end - start;


	if (wav) {
		// Pitch Detection
		cout << "Detecting Pitches ... ";
		start = high_resolution_clock::now();
		vector<double> pitches = detectPitches(dftout);
		stop = high_resolution_clock::now();
		duration = stop - start;
		cout << "done in " << duration.count() << " seconds" << endl;
		durations.push_back(duration.count());
		cout << "Pitches detected:" << endl;
		for (unsigned int i = 0; i < pitches.size(); i++)
			cout << pitches[i] << endl;
		cout << endl;
	}

	writeVectorToFile("dur.txt", durations, (int)durations.size());

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

void writeVectorToFile(string file, vector<double> vec, int N)
{
    #ifdef PRINT
		cout << "Writing to file ... ";
		ofstream f;
		f.open("out/" + file, ios::trunc);
		for (int i = 0 ; i < N; i++)
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

vector<complex<double>> DFT(vector<double> input)
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

vector<complex<double>> PDFT(vector<double> input)
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

	#pragma omp parallel for default(none) private(real, imag, innerSum, k, tmp, n) shared(input, output, K, N, test)
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
		output.at(k) = innerSum;
	}

    return output;
}

vector<complex<double>> CT(vector<double> input)
{
    vector<complex<double>> output;

	int high = ceil(log(input.size()) / log(2));
	uint32_t pad = pow(2, high) - input.size();
	vector<double>::iterator it = input.begin();
	int half = floor(input.size() / 2);
	
	input.insert(it + half, pad, (double)0.0);

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
			complex<double> t, u, tmpk, tmpm;
			int k;
			for (k = j; k < top; k += m) {
				t = (w * output[k + m2]);
				u = output[k];
				tmpk = u + t;
				tmpm = u - t;
				if (fabs(tmpk.real()) < test)
					tmpk.real(0.0);
				if (fabs(tmpk.imag()) < test)
					tmpk.imag(0.0);
				if (fabs(tmpm.real()) < test)
					tmpm.real(0.0);
				if (fabs(tmpm.imag()) < test)
					tmpm.imag(0.0);
				output[k] = tmpk;
				output[k + m2] = tmpm;
			}
			w *= wm;
		}
	}
    return output;
}

vector<complex<double>> CTP1(vector<double> input)
{
    vector<complex<double>> output;

	int high = ceil(log(input.size()) / log(2));
	uint32_t pad = pow(2, high) - input.size();
	vector<double>::iterator it = input.begin();
	int half = floor(input.size() / 2);
	
	input.insert(it + half, pad, (double)0.0);

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
			complex<double> t, u, tmpk, tmpm;
			int k;
			#pragma omp parallel for default(none) shared(j, output, test, top, m, m2, w) private(t, u, k, tmpk, tmpm)
			for (k = j; k < top; k += m) {
				t = (w * output[k + m2]);
				u = output[k];
				tmpk = u + t;
				tmpm = u - t;
				if (fabs(tmpk.real()) < test)
					tmpk.real(0.0);
				if (fabs(tmpk.imag()) < test)
					tmpk.imag(0.0);
				if (fabs(tmpm.real()) < test)
					tmpm.real(0.0);
				if (fabs(tmpm.imag()) < test)
					tmpm.imag(0.0);
				#pragma omp critical
				{
					output[k] = tmpk;
					output[k + m2] = tmpm;
				}
			}
			w *= wm;
		}
	}
    return output;
}

vector<complex<double>> CTP2(vector<double> input)
{
    vector<complex<double>> output;

	int org = input.size();
	int high = log(org) / log(2);

	if (pow(2, high) != input.size()) {
		cout << "CTP2 only takes power of 2 inputs .... sorry!" << endl;
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
			complex<double> t, u, tmpk, tmpm;
			int k;
			#pragma omp parallel for default(none) shared(j, output, test, top, m, m2, w, locks) private(t, u, k, tmpk, tmpm)
			for (k = j; k < top; k += m) {
				t = (w * output[k + m2]);
				u = output[k];
				tmpk = u + t;
				tmpm = u - t;
				if (fabs(tmpk.real()) < test)
					tmpk.real(0.0);
				if (fabs(tmpk.imag()) < test)
					tmpk.imag(0.0);
				if (fabs(tmpm.real()) < test)
					tmpm.real(0.0);
				if (fabs(tmpm.imag()) < test)
					tmpm.imag(0.0);
				omp_set_lock(&(locks[k]));
				output[k] = tmpk;
				omp_unset_lock(&(locks[k]));
				omp_set_lock(&(locks[k + m2]));
				output[k + m2] = tmpm;
				omp_unset_lock(&(locks[k + m2]));
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

