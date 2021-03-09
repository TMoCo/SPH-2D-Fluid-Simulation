#ifndef VECTOR2_H
#define VECTOR2_H

#include <iostream>

class Vector2
    { // Cartesian3
    public:
    // the coordinates
    float x, y;

    // constructors
    Vector2();
    Vector2(float X, float Y);
    Vector2(const Vector2 &other);
    
    // equality operator
    bool operator ==(const Vector2 &other) const;

    // addition operator
    Vector2 operator +(const Vector2 &other) const;

    // subtraction operator
    Vector2 operator -(const Vector2 &other) const;
    
    // multiplication operator
    Vector2 operator *(float factor) const;

    // division operator
    Vector2 operator /(float factor) const;

    // unary minus operator
    Vector2 operator -() const;

    // += operator
    Vector2& operator +=(const Vector2& other);

    // dot product routine
    float dot(const Vector2 &other) const;
    
    // routine to find the length
    float length() const;
    
    // normalisation routine
    Vector2 unit() const;

    // returns an orthogonal vector
    // counter clockwise
    Vector2 orthogonalCCW() const;
    // clockwise
    Vector2 orthogonalCW() const;
    
    // operator that allows us to use array indexing instead of variable names
    float &operator [] (const int index);
    const float &operator [] (const int index) const;
};

// multiplication operator
Vector2 operator *(float factor, const Vector2 &right);
#endif