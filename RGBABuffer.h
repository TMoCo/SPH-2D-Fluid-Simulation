#ifndef RGBABUFFER_H
#define RGBABUFFER_H

#include "RGBAValue.h"

// A rectangular buffer containing RGBAValues (float) for rendering the simulation as a texture
class RGBABuffer { 
    public:
    //  the raw data
    RGBAValue<float> *block;
    // dimensions of the buffer
    unsigned int width, height;

    // constructor
    RGBABuffer();
    // destructor
    ~RGBABuffer();
    
    // resizes the image, destroying any contents
    void Resize(long Width, long Height);

    // indexing - retrieves the beginning of a line
    // array indexing will then retrieve an element
    RGBAValue<float> * operator [](const int rowIndex);
    
    // similar routine for const pointers
    const RGBAValue<float> * operator [](const int rowIndex) const;
    
}; 
#endif