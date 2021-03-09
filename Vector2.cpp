#include "Vector2.h"
#include "math.h"
#include <stdexcept>

// constructors
Vector2::Vector2() : x(0.0), y(0.0) {}

Vector2::Vector2(float X, float Y) : x(X), y(Y) {}

Vector2::Vector2(const Vector2 &other) : x(other.x), y(other.y) {}
    
bool Vector2::operator ==(const Vector2 &other) const { 
    return ((x == other.x) && (y == other.y));
}

Vector2 Vector2::operator +(const Vector2 &other) const {
    Vector2 returnVal(x + other.x, y + other.y);
    return returnVal;
}

Vector2 Vector2::operator -(const Vector2 &other) const {
    Vector2 returnVal(x - other.x, y - other.y);
    return returnVal;
} 

Vector2 Vector2::operator *(float factor) const {
    Vector2 returnVal(x * factor, y * factor);
    return returnVal;
}

Vector2 Vector2::operator /(float factor) const {
    Vector2 returnVal(x / factor, y / factor);
    return returnVal;
}

// unary minus operator flips inverses direction
Vector2 Vector2::operator -() const {
    Vector2 returnVal(-x, -y);
    return returnVal;
}

Vector2& Vector2::operator +=(const Vector2& other) {
    this->x += other.x;
    this->y += other.y;
    return *this;
}

// dot product routine
float Vector2::dot(const Vector2 &other) const {
    float returnVal = x * other.x + y * other.y;
    return returnVal;
}

// routine to find the length
float Vector2::length() const {
    return sqrt(x*x + y*y);   
}

// normalisation routine
Vector2 Vector2::unit() const {
    float length = sqrt(x*x+y*y);
    Vector2 returnVal;
    if (length > 00.1) {
        returnVal.x = x/length;
        returnVal.y = y/length;
    }
    else {
        returnVal.x = x;
        returnVal.y = y;
    }
    return returnVal;
}

// operator that allows us to use array indexing instead of variable names
float &Vector2::operator [] (const int index){
    // use default to catch out of range indices
    switch (index) { 
        case 0:
            return x;
        case 1:
            return y;
        default:
            throw std::out_of_range ("Index out of range");
    }
}

// operator that allows us to use array indexing instead of variable names
const float &Vector2::operator [] (const int index) const {
    // use default to catch out of range indices
    switch (index) { 
        case 0:
            return x;
        case 1:
            return y;
        default:
            throw std::out_of_range ("Index out of range");
    }
}

// returns an orthogonal vector
// counterclockwise
Vector2 Vector2::orthogonalCW() const {
    Vector2 ortho(y, -x);
    return ortho;
}

// clockwise
Vector2 Vector2::orthogonalCCW() const {
    Vector2 ortho(-y, x);
    return ortho;
}

// multiplication operator
Vector2 operator *(float factor, const Vector2 &right) {
    // scalar multiplication is commutative, so flip & return
    return right * factor;
}

