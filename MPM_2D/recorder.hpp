//
//  recorder.hpp
//  MPM_2D
//
//  Created by LeeSangheon on 2024/02/19.
//

#ifndef recorder_hpp
#define recorder_hpp

#include <iostream>
#include <stdio.h>
#include <Metal/Metal.hpp>

#include "CoreFoundation/CoreFoundation.h"

#include "time.h"

using namespace std;

void save_image_from_texture(MTL::Texture* texture);

class Recorder {
public:
    Recorder(int width, int height);
    void addFrame(MTL::Texture* texture, float time);
private:
    const int width;
    const int height;
    string path;
    MTL::Region region;
    unsigned char* pixels;
};

#endif /* recorder_hpp */
