//////////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  ------------------------
//  RGBAImage.cpp
//  ------------------------
//  
//  A minimal class for an image in single-byte RGBA format
//
//  Modified by Thomas Moreno Cooper
//  
///////////////////////////////////////////////////

#include <stdlib.h>
#include <iostream>

#include "RGBABuffer.h"

// arbitrary size for a maximum image... beware large displays!
#define MAX_IMAGE_DIMENSION 4096

// constructor
RGBABuffer::RGBABuffer() : block(NULL), width(0), height(0) {}

//  destructor
RGBABuffer::~RGBABuffer() { 
    free(block);
}

// resizes the image, destroying any contents
void RGBABuffer::Resize(long Width, long Height) {
    // check validity of dimensions
    if ((Width < 0) || (Width > MAX_IMAGE_DIMENSION) || (Height < 0) || (Height > MAX_IMAGE_DIMENSION))
        throw std::runtime_error("Could not handle frame buffer size.");

    // if our old block is non-null, release the old pointer
    if (block != NULL)
        free(block);

    // use calloc() to allocate & zero memory
    block = (RGBAValue<float>*) calloc(Height * Width, sizeof (RGBAValue<float>));
    if (block == NULL)
        throw std::runtime_error("Could not allocate frame buffer memory.");

    // now that it's reallocated and copied, reset the parameters
    height = Height;
    width = Width;
}

// indexing - retrieves the beginning of a line
// array indexing will then retrieve an element
RGBAValue<float>* RGBABuffer::operator [](const int rowIndex) {
    // use pointer arithmetic to compute the row beginning
    return block+(rowIndex*width);
}
    
// similar routine for const pointers
const RGBAValue<float>* RGBABuffer::operator [](const int rowIndex) const {
    // use pointer arithmetic to compute the row beginning
    return block+(rowIndex*width);
}
