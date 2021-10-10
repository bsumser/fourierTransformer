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

    //get imaginary number i for DFT algorithm
    int I = -1;
    I = sqrt(I);

    //sigma notation algorithm start
    for (int i = 0; i < input.size() - 1; i++)
    {
        input[i] * (exp((-1 * I * 2 * M_PI) / input.size()));
    }
    return sum;
}