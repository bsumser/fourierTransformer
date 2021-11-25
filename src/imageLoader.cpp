#include "../include/imageLoader.h"
#include <iostream>
#include <string>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>


ImageLoader::ImageLoader(std::string path)
{
    std::cout << "constructor called for image loader" << std::endl;
    std::cout << "loading image with path " << path << std::endl;
}
