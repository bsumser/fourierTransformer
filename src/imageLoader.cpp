#include "../include/imageLoader.h"
#include <iostream>
#include <string>
#include "../include/lodepng.h"

ImageLoader::ImageLoader(const char* path)
{
    std::cout << "constructor called for image loader" << std::endl;
    std::cout << "loading image with path " << path << std::endl;

    std::vector<unsigned char> img;
    unsigned w, h;

    unsigned error = lodepng::decode(img, w, h, path);

    if (error) std::cout << "decoding error" << error << ": " << lodepng_error_text(error) << std::endl;
    else
    {
        image = img;
        width = w;
        height = h;
        std::cout << "Success" << std::endl;
        for (int i = 0; i < image.size(); i+=4)
        {
            std::cout << "Pixel " << i << " is : " 
            << "r:" << (int)image[i] << " " 
            << "g:" << (int)image[i + 1] << " "
            << "b:" << (int)image[i + 2] << " "
            << "a:" << (int)image[i + 3] << " "
            << std::endl;
        }
    }
}        
void ImageLoader::grayscaler()
{
    std::cout << "attempting to encode grayscale image" << std::endl;
    for (int i = 0; i < image.size(); i+=4)
    {
        char grayPixel;
        unsigned int gray = ((int)image[i] + (int)image[i + 1] + (int)image[i + 2]) / 3;
        grayPixel = gray;
        image[i] = grayPixel;
        image[i + 1] = grayPixel;
        image[i + 2] = grayPixel;
        std::cout << "Pixel " << i << " is : " 
        << "r:" << (int)image[i] << " " 
        << "g:" << (int)image[i + 1] << " "
        << "b:" << (int)image[i + 2] << " "
        << "a:" << (int)image[i + 3] << " "
        << std::endl;
    }
    char* filename = "test.png";

    unsigned error = lodepng::encode(filename, image, width, height);

    //if there's an error, display it
    if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}