//
//  recorder.cpp
//  MPM_2D
//
//  Created by LeeSangheon on 2024/02/19.
//

#include "recorder.hpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

Recorder::Recorder(int width, int height) : width(width), height(height) {

    CFBundleRef mainBundle = CFBundleGetMainBundle();
    CFURLRef resourcesURL = CFBundleCopyResourcesDirectoryURL(mainBundle);
    char path_[PATH_MAX];
    if (!CFURLGetFileSystemRepresentation(resourcesURL, TRUE, (UInt8 *)path_, PATH_MAX))
    {
        std::cout << ("Error");
    }
    CFRelease(resourcesURL);
    strcat(path_, "/VideoResult");
    path = string(path_);
    
    string command = string("mkdir -p ") + path;
    system(command.c_str());
    
    time_t rawTime = time(NULL);
    struct tm* pTimeInfo;
    
    pTimeInfo = localtime(&rawTime);
    
    string foldername = (to_string(pTimeInfo->tm_year + 1900) + "-" +
                         to_string(pTimeInfo->tm_mon + 1) + "-" +
                         to_string(pTimeInfo->tm_mday) + "-" +
                         to_string(pTimeInfo->tm_hour) + "-" +
                         to_string(pTimeInfo->tm_min) + "-" +
                         to_string(pTimeInfo->tm_sec));
    
    path += "/" + foldername;
    
    command = string("mkdir -p ") + path;
    system(command.c_str());

    //chdir(path);
    std::cout << "Frames will be saved in " << path << "\n";
        
    region = MTL::Region(0, 0, width, height);
    unsigned char _pixels[width*height*4];
    pixels = _pixels;
}

void Recorder::addFrame(MTL::Texture* texture, float time) {
    texture->getBytes(pixels, 4*width, region, 0);
    for (int i=0; i<width*height; i++) {
        char b = pixels[4*i+0];
        char g = pixels[4*i+1];
        char r = pixels[4*i+2];
        char a = pixels[4*i+3];
        pixels[4*i+0] = r;
        pixels[4*i+1] = g;
        pixels[4*i+2] = b;
        pixels[4*i+3] = a;
    }
    stbi_write_png(&((path+string("/")+to_string(time)+string(".png"))[0]),
                   width, height, 4, &(pixels[0]), 4*width);
}
