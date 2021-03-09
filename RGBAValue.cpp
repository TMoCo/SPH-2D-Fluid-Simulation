#include "RGBAValue.h"

template class RGBAValue<int>;
template class RGBAValue<float>;
template class RGBAValue<double>;

// default constructor
template <typename T>
RGBAValue<T>::RGBAValue() : red(0), green(0), blue(0), alpha(0) {}

// value constructor with floats
// values outside 0.0-255.0 are clamped
template <typename T>
RGBAValue<T>::RGBAValue(T Red, T Green, T Blue, T Alpha) { 
    if (Red > 255.0) Red = 255.0;
    if (Red < 0.0) Red = 0.0;
    // assign and let the compiler convert to unsigned char
    red = Red;

    if (Green > 255.0) Green = 255.0;
    if (Green < 0.0) Green = 0.0;
    green = Green;

    if (Blue > 255.0) Blue = 255.0;
    if (Blue < 0.0) Blue = 0.0;
    blue = Blue;

    if (Alpha > 255.0) Alpha = 255.0;
    if (Alpha < 0.0) Alpha = 0.0;
    alpha = Alpha;
} 

// copy constructor
template <typename T>
RGBAValue<T>::RGBAValue(const RGBAValue<T>& other) :
    red(other.red), green(other.green), blue(other.blue), alpha(other.alpha) {}

template <typename T>
RGBAValue<T> RGBAValue<T>::modulate(const RGBAValue<T> &right) const {
    T leftRed = (T) red / 255.0, leftGreen = (T) green / 255.0, 
        leftBlue = (T) blue / 255.0;
    T rightRed = (T) right.red / 255.0, rightGreen = (T) right.green / 255.0, 
        rightBlue = (T) right.blue / 255.0;
    return RGBAValue<T>(255 * leftRed * rightRed, 255 * leftGreen * rightGreen, 
        255 * leftBlue * rightBlue, 255);
}

template <typename T>
RGBAValue<T>& RGBAValue<T>::operator=(const RGBAValue<T>& copy) {
    red = copy.red;
    green = copy.green;
    blue = copy.blue;
    alpha = copy.alpha;
    return *this;
}

