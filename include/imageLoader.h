#ifndef IMAGELOADER_H
#define IMAGELOADER_H
#include <string>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>


class ImageLoader
{
    public:
        ImageLoader(std::string path);
    private:
        std::string path;
};

#endif