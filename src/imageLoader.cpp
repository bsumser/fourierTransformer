#include "../include/imageLoader.h"
#include <iostream>
#include <string>
#include <lodepng.h>


ImageLoader::ImageLoader(std::string path)
{
    std::cout << "constructor called for image loader" << std::endl;
    std::cout << "loading image with path " << path << std::endl;
}
