#ifndef IMAGELOADER_H
#define IMAGELOADER_H
#include <string>
#include "../include/lodepng.h"


class ImageLoader
{
    public:
        ImageLoader(const char* path);
        void grayscaler();
        const char* path;
        std::vector<unsigned char> image;
        unsigned width, height;
    private:
};

#endif