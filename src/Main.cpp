#include <iostream>
#include <vector>
#include <cmath>
#include "../include/Timelord.h"

int main() {
    std::cout << "Hello world" << std::endl;

    Timelord newTimeLord();
    return 0;
}

int discreteFourierTransform(std::vector<float> input)
{
    int sum = 0;
    int inputSize = input.size();

    //get imaginary number i for DFT algorithm
    int eulerPower = -1;
    eulerPower = sqrt(eulerPower) * 2 * M_PI;
    eulerPower = (eulerPower / inputSize) * -1;

    //sigma notation algorithm start
    for (int i = 0; i < input.size() - 1; i++)
    {
        input[i] * exp(eulerPower);
    }
    return sum;
}