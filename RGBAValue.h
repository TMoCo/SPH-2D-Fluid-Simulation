#ifndef RGBAVALUE_H
#define RGBAVALUE_H

template <typename T>
class RGBAValue {
    public:
    // just a container for the components
    unsigned char red, green, blue, alpha;

    // default constructor
    RGBAValue();

    // value constructor with template T
    // values outside 0.0-255.0 are clamped
    RGBAValue(T Red, T Green, T Blue, T Alpha);

    // copy constructor
    RGBAValue(const RGBAValue& other);

    RGBAValue modulate(const RGBAValue& right) const;

    // assignment operator
    RGBAValue& operator =(const RGBAValue& copy);
    
};

#endif
