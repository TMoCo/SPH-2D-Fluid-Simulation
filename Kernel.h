#ifndef KERNEL_H
#define KERNEL_H

#include "Vector2.h"

// pi for computing kernel values
#define PI 3.14159265358979323846

// some small threshold value
#define EPSILON 0.001

// the definition of the kernels in a static Kernel class
class Kernel {
    public:
    // Poly6 smoothing kernel
    static float Poly6(const Vector2& r, const float& h);
    // Poly6 gradient
    static float Poly6Gradient(const Vector2& r, const float& h);
    // Poly6 laplacian
    static float Poly6Laplacian(const Vector2& r, const float& h);

    // Spiky smoothing kernel
    static float Spiky(const Vector2& r, const float& h);
    // Spiky gradient
    static float SpikyGradient(const Vector2& r, const float& h);

    // Viscosity smoothing kernel
    static float Viscosity(const Vector2& r, const float& h);
    // Viscosity gradient
    static float ViscosityGradient(const Vector2& r, const float& h);
    // Viscosity laplacian
    static float ViscosityLaplacian(const Vector2& r, const float& h);

    // purely static class
    private:
    Kernel();

};

#endif