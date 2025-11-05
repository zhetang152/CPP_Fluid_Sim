#ifndef MATH_H
#define MATH_H

#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <string>

#ifdef FLOAT_AS_Float
using Float = double;
#else
using Float = float;
#endif

constexpr Float Pi = 3.14159265358979323846;

constexpr Float PosInfinity = std::numeric_limits<Float>::infinity();
constexpr Float NegInfinity = -std::numeric_limits<Float>::infinity();

class Interval{
private:
    Float low;
    Float high;
public:
    Interval() = default;
    explicit Interval(Float l, Float h) : low(std::min(l,h)), high(std::max(l,h)) {}
    Interval(Float c) : low(c), high(c) {}

    Float LowerBound() const { return low; }
    Float UpperBound() const { return high; }
    Float Midpoint() const { return (low + high) / 2; }
    Float Width() const { return high - low; }

    Interval operator+(const Interval& other) const {
        return Interval(AddRoundDown(low, other.low), AddRoundUp(high, other.high));
    }
    Interval operator-(const Interval& other) const {
        return Interval(SubRoundDown(low, other.high), SubRoundUp(high, other.low));
    }
    Interval operator*(const Interval& other) const {
        Float down[4]= {MulRoundDown(low, other.low), MulRoundDown(low, other.high),
                         MulRoundDown(high, other.low),MulRoundDown(high, other.high)};
        Float up[4]= {MulRoundUp(low, other.low), MulRoundUp(low, other.high),
                       MulRoundUp(high, other.low),MulRoundUp(high, other.high)};
        return Interval(*std::min_element(down, down+4), *std::max_element(up, up+4));
    }
};

inline Float AddRoundUp(Float a, Float b) {
    return std::nextafter(a + b, PosInfinity);
}
inline Float AddRoundDown(Float a, Float b) {
    return std::nextafter(a+b, NegInfinity);
}
inline Float SubRoundUp(Float a, Float b) {
    return std::nextafter(a - b, PosInfinity);
}
inline Float SubRoundDown(Float a, Float b) {
    return std::nextafter(a - b, NegInfinity);
}
inline Float MulRoundUp(Float a, Float b) {
    return std::nextafter(a * b, PosInfinity);
}
inline Float MulRoundDown(Float a, Float b) {
    return std::nextafter(a * b, NegInfinity);
}
inline Float DivRoundUp(Float a, Float b) {
    return std::nextafter(a / b, PosInfinity);
}
inline Float DivRoundDown(Float a, Float b) {
    return std::nextafter(a / b, NegInfinity);
}


template <typename T>
inline constexpr T Sqr(T v) {
    return v * v;
}


inline Float SafeAsin(Float x){

}

template <typename Ta, typename Tb, typename Tc, typename Td>
inline auto DifferenceOfProducts(Ta a, Tb b, Tc c, Td d) {
    auto cd = c * d;
    auto differenceOfProducts = FMA(a, b, -cd);
    auto error = FMA(-c, d, cd);
    return differenceOfProducts + error;
}

#endif // MATH_H