#ifndef IMAGELOADER_H
#define IMAGELOADER_H
#include <string>
#include <lodepng.h>


class ImageLoader
{
    public:
        ImageLoader(std::string path);
    private:
        std::string path;
};

#endif