#ifndef IMAGELOADER_H
#define IMAGELOADER_H
#include <string>
#include "../include/lodepng.h"


class ImageLoader
{
    public:
        ImageLoader(const char* path);
    private:
        std::string path;
};

#endif