#include "Kernel.h"
#include <math.h>


float Kernel::Poly6(const Vector2& r, const float& h) {
    float distance = r.length();
    // return 0 if the length of r is outside of the range
    if ((distance < EPSILON) || (distance > h))
        return 0;
    return 315.0 / (64.0 * PI * pow(h, 9)) * pow(pow(h, 2) - pow(distance, 2), 3);
}

float Kernel::Poly6Gradient(const Vector2& r, const float& h) {
    float distance = r.length();
    // return 0 if the length of r is outside of the range
    if ((distance < EPSILON) || (distance > h))
        return 0;
    return -945.0 * distance * pow((pow(h,2) - pow(distance,2)), 2) / (32.0 * PI * pow(h,9)); 
}

float Kernel::Poly6Laplacian(const Vector2& r, const float& h) {
    float distance = r.length();
    // return 0 if the length of r is outside of the range
    if ((distance < EPSILON) || (distance > h))
        return 0;
    return -945.0 * (pow(h,2) - 5*pow(distance,2)) * (pow(h,2) - pow(distance,2)) / (32.0 * PI * pow(h,9));
}

float Kernel::Spiky(const Vector2& r, const float& h) {
    float distance = r.length();
    // return 0 if the length of r is outside of the range
    if ((distance < EPSILON) || (distance > h))
        return 0;
    return 15.0 / (PI * pow(h, 6)) * pow(h - distance, 3);
}

float Kernel::SpikyGradient(const Vector2& r, const float& h) {
    float distance = r.length();
    // return 0 if the length of r is outside of the range
    if ((distance < EPSILON) || (distance > h))
        return 0;
    return -45.0 * pow(h - distance, 2) / (PI * pow(h, 6));
}

float Kernel::Viscosity(const Vector2& r, const float& h) {
    float distance = r.length();
    // return 0 if the length of r is outside of the range
    if ((distance < EPSILON) || (distance > h))
        return 0;
    return -15.0 * (-pow(distance,3) / (2*pow(h,2)) + pow(distance,2) / 
        pow(h,2) + h / (2.0*distance) - 1) / (2.0 * PI * pow(h, 3));
}

float Kernel::ViscosityLaplacian(const Vector2& r, const float& h) {
    float distance = r.length();
    // return 0 if the length of r is outside of the range
    if ((distance < EPSILON) || (distance > h))
        return 0;
    return 45.0 / (PI * pow(h,6)) * (h - distance);
}



