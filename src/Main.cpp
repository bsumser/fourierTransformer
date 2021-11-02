#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include "../include/Timelord.h"

using namespace std::chrono;


std::vector<float> floatGenerator(int vectorSize);
int discreteFourierTransform(std::vector<float> input);

int main(int argc, char* argv[]) {
    Timelord newTimeLord();
    floatGenerator(100000);
    return 0;
}

int discreteFourierTransform(std::vector<float> input)
{
    float sum = 0;
    int inputSize = input.size();

    //get imaginary number i for DFT algorithm
    float eulerPower = (sqrt(-1) * 2 * M_PI) / inputSize;

    //sigma notation algorithm start
    for (int i = 0; i < input.size() - 1; i++)
    {
        input[i] * exp(eulerPower);
    }
    return sum;
}

std::vector<float> floatGenerator(int vectorSize)
{
    auto start = high_resolution_clock::now();
    std::default_random_engine floatGenerator;
    std::mt19937 mt(floatGenerator());
    std::uniform_real_distribution<float> uniform_distance(0, 10.001);
    std::vector<float> floatVector;

    for (int i = 0; i < vectorSize - 1; i++) {
        floatVector.push_back(uniform_distance(mt));
        std::cout << "rand is " << uniform_distance(mt) << std::endl;
    }
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    return floatVector;
}