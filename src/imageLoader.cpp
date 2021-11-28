#include "../include/imageLoader.h"
#include <iostream>
#include <string>
#include "../include/lodepng.h"

std::vector<unsigned char> img;
unsigned w, h;

ImageLoader::ImageLoader(const char* path)
{
    std::cout << "constructor called for image loader" << std::endl;
    std::cout << "loading image with path " << path << std::endl;

    std::vector<unsigned char> image;
    unsigned width, height;

    unsigned error = lodepng::decode(image, width, height, path);

    if (error) std::cout << "decoding error" << error << ": " << lodepng_error_text(error) << std::endl;
    else
    {
        img = image;
        w = width;
        h = height;
        std::cout << "Success" << std::endl;
        for (int i = 0; i < img.size(); i+=4)
        {
            std::cout << "Pixel " << i << " is : " 
            << "r:" << (int)img[i] << " " 
            << "g:" << (int)img[i + 1] << " "
            << "b:" << (int)img[i + 2] << " "
            << "a:" << (int)img[i + 3] << " "
            << std::endl;
        }

        std::cout << "attempting to encode grayscale image" << std::endl;
        for (int i = 0; i < img.size(); i+=4)
        {
            char grayPixel;
            unsigned int gray = ((int)img[i] + (int)img[i + 1] + (int)img[i + 2]) / 3;
            grayPixel = gray;
            img[i] = grayPixel;
            img[i + 1] = grayPixel;
            img[i + 2] = grayPixel;
            std::cout << "Pixel " << i << " is : " 
            << "r:" << (int)img[i] << " " 
            << "g:" << (int)img[i + 1] << " "
            << "b:" << (int)img[i + 2] << " "
            << "a:" << (int)img[i + 3] << " "
            << std::endl;
        }
        char* filename = "test.png";

        unsigned error = lodepng::encode(filename, img, w, h);

        //if there's an error, display it
        if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
    }
}
