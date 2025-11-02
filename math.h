#pragma once
#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <string>
#include "pbrt.h"


constexpr Float Pi = 3.14159265358979323846;

constexpr double PosInfinity = std::numeric_limits<double>::infinity();
constexpr double NegInfinity = -std::numeric_limits<double>::infinity();

class Interval{
private:
    double low;
    double high;
public:
    Interval() = default;
    explicit Interval(double l, double h) : low(std::min(l,h)), high(std::max(l,h)) {}
    Interval(double c) : low(c), high(c) {}

    double LowerBound() const { return low; }
    double UpperBound() const { return high; }
    double Midpoint() const { return (low + high) / 2; }
    double Width() const { return high - low; }

    Interval operator+(const Interval& other) const {
        return Interval(AddRoundDown(low, other.low), AddRoundUp(high, other.high));
    }
    Interval operator-(const Interval& other) const {
        return Interval(SubRoundDown(low, other.high), SubRoundUp(high, other.low));
    }
    Interval operator*(const Interval& other) const {
        double down[4]= {MulRoundDown(low, other.low), MulRoundDown(low, other.high),
                         MulRoundDown(high, other.low),MulRoundDown(high, other.high)};
        double up[4]= {MulRoundUp(low, other.low), MulRoundUp(low, other.high),
                       MulRoundUp(high, other.low),MulRoundUp(high, other.high)};
        return Interval(*std::min_element(down, down+4), *std::max_element(up, up+4));
    }
};

inline double AddRoundUp(double a, double b) {
    return std::nextafter(a + b, PosInfinity);
}
inline double AddRoundDown(double a, double b) {
    return std::nextafter(a+b, NegInfinity);
}
inline double SubRoundUp(double a, double b) {
    return std::nextafter(a - b, PosInfinity);
}
inline double SubRoundDown(double a, double b) {
    return std::nextafter(a - b, NegInfinity);
}
inline double MulRoundUp(double a, double b) {
    return std::nextafter(a * b, PosInfinity);
}
inline double MulRoundDown(double a, double b) {
    return std::nextafter(a * b, NegInfinity);
}
inline double DivRoundUp(double a, double b) {
    return std::nextafter(a / b, PosInfinity);
}
inline double DivRoundDown(double a, double b) {
    return std::nextafter(a / b, NegInfinity);
}


template <typename T>
inline constexpr T Sqr(T v) {
    return v * v;
}


inline float SafeAsin(float x){

}

template <typename Ta, typename Tb, typename Tc, typename Td>
inline auto DifferenceOfProducts(Ta a, Tb b, Tc c, Td d) {
    auto cd = c * d;
    auto differenceOfProducts = FMA(a, b, -cd);
    auto error = FMA(-c, d, cd);
    return differenceOfProducts + error;
}