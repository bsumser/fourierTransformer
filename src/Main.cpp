#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <chrono>
#include "../include/Timelord.h"

using namespace std::chrono;
using namespace std;

vector<float> floatGenerator(int vectorSize);
void processArgs(int argc, char* argv[]);
vector<std::complex<double>> discreteFourierTransform(vector<complex<double>> input);
vector<std::complex<double>> discreteFourierTransformFaster(vector<std::complex<double>> input);
vector<std::complex<double>> discreteFourierTransformTurkey(vector<std::complex<double>> input);
vector<std::complex<double>> signalGenerator(int sampleSize);

int main(int argc, char* argv[]) {
    processArgs(argc, argv);
    Timelord newTimeLord();
    //floatGenerator(100000);

    vector<std::complex<double>> input = signalGenerator(1000);
    discreteFourierTransform(input);
    discreteFourierTransformFaster(input);
    discreteFourierTransformTurkey(input);

    return 0;
}

void processArgs(int argc, char* argv[])
{
    cout << "processing args" << endl;
    for (int i = 0; i < argc; i++) {
        cout << "arg " << i << " is " << argv[i] << endl;
    }
}

vector<std::complex<double>> discreteFourierTransform(vector<std::complex<double>> input)
{
    //initialize sizes of samples
    int N = input.size();
    int K = input.size();

    //init variable for internal loop
    std::complex<double> innerSum;

    //init vector of std::complex doubles for results
    vector<std::complex<double>> output;

    //set output vector to have enough space
    output.reserve(N);

    auto start = high_resolution_clock::now();  //start clock to measure execution time
    //sigma notation algorithm start
    for (int k = 0; k < K; k++)
    {
        innerSum = std::complex<double>(0,0);
        for (int n = 0; n < N; n++) {
            //process real and imaginary parts of sum via definition in "https://en.wikipedia.org/wiki/Discrete_Fourier_transform"
            double real = cos(((2 * M_PI) / N) * k * n);
            double imag = sin(((2 * M_PI) / N) * k * n);

            //store as std::complex double for multiplication
            std::complex<double> w (real, -imag);
            innerSum += input[n] * w;
        }
        //add value to the output vector
        output.push_back(innerSum);
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << __func__ << " returned in " << duration.count() << " microseconds" << endl;

    for (int i = 0; i < 6; i++) {
        cout << output[i] << endl;
    }
    
    return output;
}

vector<std::complex<double>> discreteFourierTransformFaster(vector<std::complex<double>> input)
{
    //initiliaze sizes of samples
    int N = input.size();
    int K = input.size();

    //loop indexes
    int k;
    int n;

    //init variable for internal loop
    std::complex<double> innerSum;

    //init vector of std::complex doubles for results
    vector<std::complex<double>> output;

    //set output vector to have enough space
    output.reserve(N);

    auto start = high_resolution_clock::now();  //start clock to measure execution time

    for (k = 0; k < K; k++)
    {
        innerSum = std::complex<double>(0,0);
        //TODO fix this reduction
        for (n = 0; n < N; n++) {
            //process real and imaginary parts of sum via definition in "https://en.wikipedia.org/wiki/Discrete_Fourier_transform"
            double real = cos(((2 * M_PI) / N) * k * n);
            double imag = sin(((2 * M_PI) / N) * k * n);

            //store as std::complex double for multiplication
            std::complex<double> w (real, -imag);
            innerSum += input[n] * w;
        }
        //add value to the output vector
        output.push_back(innerSum);
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << __func__ << " returned in " << duration.count() << " microseconds" << endl;

    for (int i = 0; i < 6; i++) {
        cout << output[i] << endl;
    }
    
    return output;
}

vector<std::complex<double>> discreteFourierTransformTurkey(vector<std::complex<double>> input)
{
    //initiliaze sizes of samples
    int N = input.size();
    int K = input.size();

    //init variable for internal loop
    std::complex<double> innerSum;

    //init vector of std::complex doubles for results
    vector<std::complex<double>> output;

    //set output vector to have enough space
    output.reserve(N);

    auto start = high_resolution_clock::now();  //start clock to measure execution time
    
    for (int k = 0; k < K; k++)
    {
        innerSum = std::complex<double>(0,0);
        for (int n = 0; n < N / 2 - 1; n++) {
            //process real and imaginary parts of sum via definition in "https://en.wikipedia.org/wiki/Discrete_Fourier_transform"
            double real = cos(((2 * M_PI) / (N / 2)) * k * (2 * n));
            double imag = sin(((2 * M_PI) / (N / 2)) * k * (2 * n));

            //store as std::complex double for multiplication
            std::complex<double> w (real, -imag);
            innerSum += input[n] * w;
        }
        for (int n = 1; n < N / 2 - 1; n++) {
            //process real and imaginary parts of sum via definition in "https://en.wikipedia.org/wiki/Discrete_Fourier_transform"
            double real = cos(((2 * M_PI) / (N / 2)) * k * (2 * n + 1));
            double imag = sin(((2 * M_PI) / (N / 2)) * k * (2 * n + 1));

            //store as std::complex double for multiplication
            std::complex<double> w (real, -imag);
            innerSum += input[n] * w;
        }
        //add value to the output vector
        output.push_back(innerSum);
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << __func__ << " returned in " << duration.count() << " microseconds" << endl;

    for (int i = 0; i < 6; i++) {
        cout << output[i] << endl;
    }
    
    return output;
}

vector<float> floatGenerator(int vectorSize)
{
    auto start = high_resolution_clock::now();  //start clock to measure execution time

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
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "generated vector of " << vectorSize << " elements in " << duration.count() << " microseconds" << endl;
    return floatVector;
}

vector<std::complex<double>> signalGenerator(int sampleSize)
{
    int N = sampleSize;
    double amplitude = 3;  //amplitude, peak deviation from 0
    double freq = 0;   //ordinary frequency, oscillations
    double time = 0;   //time
    double shift = M_PI / 2;  //phase shift of signal

    vector<std::complex<double>> output; //output vector to store the test signal
    output.reserve(N);  //allocate proper size for vector

    auto start = high_resolution_clock::now();  //start clock to measure execution time
    for (int i = 0; i < N; i++) {
        auto sample = std::complex<double>(cos((2 * M_PI/ N) * amplitude * i + shift),0.0);

        output.push_back(sample);
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << __func__ << " returned in " << duration.count() << " microseconds" << endl;
    return output;
}