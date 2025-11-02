#pragma once

#include <cmath>
#include <stdexcept>

#include "math.h"
#include "pbrt.h"
//类型萃取
template<typename T>
struct TLength{
    using type = Float;
};
template<>
struct TLength<double>{
    using type = double;
};
template<>
struct TLength<long double>{
    using type = long double;
};
template<>
struct TLength<Interval>{
    using type = Interval;
};

template<template<typename> class Child, typename T>
class Tuple3 {
public:
    T x{}, y{}, z{};

    Tuple3() = default;
    Tuple3(T x_val, T y_val, T z_val): x(x_val), y(y_val), z(z_val) {}
    Tuple3(const Tuple3& other): x(other.x), y(other.y), z(other.z) {}

    bool HasNaN()const{
        return std::isnan(x) || std::isnan(y) || std::isnan(z);
    }

    Tuple3 operator[](int i)const{
        if(i == 0) return x;
        if(i == 1) return y;
        if(i == 2) return z;
        throw std::out_of_range("index out of range.");
    }
    //linear operations and boolean operations
    template<typename U>
    auto operator+(Child<U> c) const -> Child<decltype(T{} + U{})>{
        return {x+c.x, y+c.y, z+c.z};
    }

    template<typename U>
    Child<T>& operator+=(const Child<U>& other){
        x += other.x;
        y += other.y;
        z += other.z;
        return static_cast<Child<T>&>(*this);
    }

    template<typename U>
    auto operator-(Child<U>& c) const -> Child<decltype(T{} - U{})>{
        return {x - c.x, y - c.y, z - c.z};
    }

    template<typename U>
    Child<T>& operator-=(const Child<U>& other){
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return static_cast<Child<T>&>(*this);
    }

    template<typename U>
    auto operator*(U s) const -> Child<decltype(T{} * U{})>{
        return {x * s, y * s, z * s};
    }

    template<typename U>
    Child<T>& operator*=(U scalar){
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return static_cast<Child<T>&>(*this);
    }

    template <typename U>
    auto operator/(U s) const -> Child<decltype(T{} / U{})>{
        return {x / s, y / s, z / s};
    }

    template<typename U>
    Child<T>& operator/=(U scalar){
        x /= scalar;
        y /= scalar;
        z /= scalar;
        return static_cast<Child<T>&>(*this);
    }

    Child<T> operator-() const { return {-x, -y, -z}; }

    bool operator==(Child<T> c) const { return x == c.x && y == c.y && z == c.z; }

    bool operator!=(Child<T> c) const {return x != c.x || y != c.y || z != c.z; }
};
//general algebraic operations
template<template<typename> class Child, typename T, typename U>
inline auto operator*(U s, Tuple3<Child, T>& t) -> Child<decltype(T{} * U{})>{
    return t * s;
}
template<template<typename> class Child, typename T>
inline Child<T> Abs(Tuple3<Child, T>& t){
    using std::abs;
    return {abs(t.x), abs(t.y), abs(t.z)};
}
template<template<typename> class Child, typename T>
inline Child<T> Ceil(Tuple3<Child, T>& t){
    using std::ceil;
    return {ceil(t.x), ceil(t.y), ceil(t.z)};
}
template<template<typename> class Child, typename T>
inline Child<T> Floor(Tuple3<Child, T>& t){
    using std::floor;
    return {floor(t.x), floor(t.y), floor(t.z)};
}
template<template<typename> class Child, typename T, typename U>
inline auto LI(U t, Tuple3<Child, T>& other1, Tuple3<Child, T>& other2){
    return (1 - t) * other1 + t * other2;
}
template<template<typename> class Child, typename T>
inline Child<T> FMA(T t, Tuple3<Child, T>&other1, Tuple3<Child, T>&other2){
    using std::fma;
    return {fma(t, other1.x, other2.x),
            fma(t, other1.y, other2.y),
            fma(t, other1.z, other2.z)};
}
template <template <class> class Child, typename T>
inline Child<T> FMA(Tuple3<Child, T> other1, T t, Tuple3<Child, T> other2) {
    return FMA(t, other1, other2);
}
//variants analysis
template <template <class> class Child, typename T>
inline int MinComponentIndex(Tuple3<Child, T> t) {
    return (t.x < t.y) ? ((t.x < t.z) ? 0 : 2) : ((t.y < t.z) ? 1 : 2);
}
template <template <class> class Child, typename T>
inline int MaxComponentIndex(Tuple3<Child, T> t) {
    return (t.x > t.y) ? ((t.x > t.z) ? 0 : 2) : ((t.y > t.z) ? 1 : 2);
}
template <template <class> class Child, typename T>
inline T MinComponentValue(Tuple3<Child, T> t) {
    using std::min;
    return min({t.x, t.y, t.z});
}
template <template <class> class Child, typename T>
inline T MaxComponentValue(Tuple3<Child, T> t) {
    using std::max;
    return max({t.x, t.y, t.z});
}
template <template <class> class Child, typename T>
inline Child<T> Min(Tuple3<Child, T> a, Tuple3<Child, T> b) {
    using std::min;
    return {min(a.x, b.x), min(a.y, b.y), min(a.z, b.z)};
}
template <template <class> class Child, typename T>
inline Child<T> Max(Tuple3<Child, T> a, Tuple3<Child, T> b) {
    using std::max;
    return {max(a.x, b.x), max(a.y, b.y), max(a.z, b.z)};
}
//unusual operations
template <template <class> class Child, typename T>
inline Child<T> Permute(Tuple3<Child, T> t, std::array<int, 3> p) {
    return {t[p[0]], t[p[1]], t[p[2]]};
}
template <template <class> class C, typename T>
inline T HProd(Tuple3<C, T> t) {
    return t.x * t.y * t.z;
}

template<typename T>
class Vector3: public Tuple3<Vector3, T> {
public:
    using Tuple3<Vector3, T>::x;
    using Tuple3<Vector3, T>::y;
    using Tuple3<Vector3, T>::z;

    Vector3() = default;
    Vector3(T x_val, T y_val, T z_val): Tuple3<Vector3, T>(x_val, y_val, z_val) {}

    template<typename U>
    explicit Vector3(Vector3<U> other): Tuple3<Vector3, T>(T(other.x), T(other.y), T(other.z)) {}
    template<typename U>
    explicit Vector3(Point3<U> p) {}
    template<typename U>
    explicit Vector3(Normal3<U> n) {}
};

using Vector3f = Vector3<Float>;
using Vector3i = Vector3<int>;

template<typename T>
class Point3: public Tuple3<Point3, T>{
public:
    using Tuple3<Point3, T>::x;
    using Tuple3<Point3, T>::y;
    using Tuple3<Point3, T>::z;

    Point3() = default;
    Point3(T x_val, T y_val, T z_val): Tuple3<Point3, T>(x_val, y_val, z_val) {}
    template<typename U>
    explicit Point3(Point3<U> other): Tuple3<Point3, T>(T(other.x), T(other.y), T(other.z)) {}
    template<typename U>
    explicit Point3(Point3<U> other): Tuple3<Vector3, T>(T(other.x), T(other.y), T(other.z)) {}

    using Tuple3<Point3, T>::operator*;
    using Tuple3<Point3, T>::operator*=;

    Point3<T> operator-() const { return {-x, -y, -z}; }

    template<typename U>
    auto operator+(Vector3<U> v)const -> Point3<decltype(T{} + U{})>{
        return {x + v.x, y + v.y, z + v.z};
    }
    template<typename U>
    Point3<T>& operator+=(Vector3<U> v){
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    template<typename U>
    auto operator-(Vector3<U> v) -> Point3<decltype(T{} - U{})>{
        return {x - v.x, y - v.y, z - v.z};
    }
    template<typename U>
    auto operator-(Point3<U> p) -> Vector3<decltype(T{} - U{})>{
        return {x - p.x, y - p.y, z - p.z};
    }
    template<typename U>
    Point3<T>& operator-=(Vector3<U> v){
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
};
template <typename T>
class Normal3 : public Tuple3<Normal3, T> {
  public:
    using Tuple3<Normal3, T>::x;
    using Tuple3<Normal3, T>::y;
    using Tuple3<Normal3, T>::z;

    using Tuple3<Normal3, T>::HasNaN;
    using Tuple3<Normal3, T>::operator+;
    using Tuple3<Normal3, T>::operator*;
    using Tuple3<Normal3, T>::operator*=;

    Normal3() = default;
    Normal3(T x, T y, T z) : Tuple3<Normal3, T>(x, y, z) {}
    template <typename U>
    explicit Normal3<T>(Normal3<U> v): Tuple3<Normal3, T>(T(v.x), T(v.y), T(v.z)) {}

    template <typename U>
    explicit Normal3<T>(Vector3<U> v): Tuple3<Normal3, T>(T(v.x), T(v.y), T(v.z)) {}
};

//operations on Vector3
//length and normalization
template<typename T>
inline T LengthSquared(const Vector3<T> v){
    return Sqr(v.x) + Sqr(v.y) + Sqr(v.z);
}
template<typename T>
inline auto Length(const Vector3<T> v) -> typename TLength<T>::type{//?
    using std::sqrt;
    return sqrt(LengthSquared(v));
}
template<typename T>
inline auto Normalize(const Vector3<T> v) -> Vector3<typename TLength<T>::type>{
    return v / Length(v);
}
//dot, angle, cross
template<typename T>
inline T Dot(const Vector3<T> v1, const Vector3<T> v2){
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
template<typename T>
inline T AbsDot(const Vector3<T> v1, const Vector3<T> v2){
    return std::abs(Dot(v1, v2));
}
template<typename T>
inline Float AngleBetween(const Vector3<T> v1, const Vector3<T> v2){
    if(Dot(v1,v2)<=0){
        return Pi - 2 * std::acos(Length(v1 + v2) / 2);
    }
    else{
        return 2 * std::acos(Length(v1 - v2) / 2);
    }
}
template<typename T>
inline Vector3<T> Cross(const Vector3<T> v1, const Vector3<T> v2){
    return {DifferenceOfProducts(v1.y, v2.z, v1.z, v2.y),//实现看一下
            DifferenceOfProducts(v1.z, v2.x, v1.x, v2.z),
            DifferenceOfProducts(v1.x, v2.y, v1.y, v2.x)};
}
template<typename T>
inline Vector3<T> GramSchmidt(const Vector3<T>& v1, const Vector3<T> v2){
    return v1 - Dot(v1, v2) * v2;
}
template<typename T>
inline void CoordinateSystem(const Vector3<T> v1, const Vector3<T> *v2, const Vector3<T> *v3){//?为什么他们用的没有const?
    Float sign = std::copysign(Float(1), v1.z);//实现看一下
    Float a = -1 / (sign + v1.z);
    Float b = v1.x * v1.y * a;
    *v2 = Vector3f(1 + sign * Sqr(v1.x) * a, sign * b, -sign * v1.x);
    *v3 = Vector3f(b, sign + Sqr(v1.y) * a, -v1.y);
}

//Operations on Point3
template<typename T>
inline auto DistanceSquared(const Point3<T> p1, const Point3<T> p2){
    return LengthSquared(p1 - p2)
}
template<typename T>
inline auto Distance(const Point3<T> p1, const Point3<T> p2){
    return Length(p1 - p2);
}

//Operations on Normal3
template<typename T>
inline auto LengthSquared(Normal3<T> n) -> typename TLength<T>::type {
    return Sqr(n.x) + Sqr(n.y) + Sqr(n.z);
}

template<typename T>
inline auto Length(Normal3<T> n) -> typename TLength<T>::type {
    using std::sqrt;
    return sqrt(LengthSquared(n));
}
template<typename T>
inline auto Normalize(Normal3<T> n) {
    return n / Length(n);
}
template <typename T>
inline auto Dot(Normal3<T> n, Vector3<T> v) -> typename TLength<T>::type {
    return n.x * v.x + n.y * v.y + n.z * v.z;
}
template <typename T>
inline auto Dot(Vector3<T> v, Normal3<T> n) -> typename TLength<T>::type {
    return Dot(n, v);
}
template <typename T>
inline auto Dot(Normal3<T> n1, Normal3<T> n2) -> typename TLength<T>::type {
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}
template <typename T>
inline auto AbsDot(Normal3<T> n, Vector3<T> v) -> typename TLength<T>::type {
    return std::abs(Dot(n, v));
}
template <typename T>
inline auto AbsDot(Vector3<T> v , Normal3<T> n) -> typename TLength<T>::type {
    return std::abs(Dot(n, v));
}
template <typename T>
inline auto AbsDot(Normal3<T> n1 , Normal3<T> n2) -> typename TLength<T>::type {
    return std::abs(Dot(n1, n2));
}
template <typename T>
inline Normal3<T> FaceForward(Normal3<T> n, Vector3<T> v) {
    return (Dot(n, v) < 0.f) ? -n : n;//?
}

template <typename T>
inline Normal3<T> FaceForward(Normal3<T> n, Normal3<T> n2) {
    return (Dot(n, n2) < 0.f) ? -n : n;
}

template <typename T>
inline Vector3<T> FaceForward(Vector3<T> v, Vector3<T> v2) {
    return (Dot(v, v2) < 0.f) ? -v : v;
}

template <typename T>
inline Vector3<T> FaceForward(Vector3<T> v, Normal3<T> n2) {
    return (Dot(v, n2) < 0.f) ? -v : v;
}

template<typename T>
class Bounds3 {
public:
    Point3<T> pMin, pMax;

    Bounds3() {
        T minNum = std::numeric_limits<T>::lowest();
        T maxNum = std::numeric_limits<T>::max();
        pMin = Point3<T>(maxNum, maxNum, maxNum);
        pMax = Point3<T>(minNum, minNum, minNum);
    }

    explicit Bounds3(const Point3<T>& p) : pMin(p), pMax(p) {}

    Bounds3(const Point3<T>& p1, const Point3<T>& p2)
        : pMin{Min(p1, p2)}, pMax{Max(p1, p2)} {}

    Point3<T> operator[](int i) const {
        if (i == 0) return pMin;
        if (i == 1) return pMax;
        throw std::out_of_range("index out of range.");
    }
    Point3<T> &operator[](int i) {
        if (i == 0) return pMin;
        if (i == 1) return pMax;
        throw std::out_of_range("index out of range.");
    }

    //坐标与点辅助函数
    Point3<T> Corner(int corner) const {
        if(corner == 0) return Point3<T>(pMin.x, pMin.y, pMin.z);
        else if(corner == 1) return Point3<T>(pMax.x, pMin.y, pMin.z);
        else if(corner == 2) return Point3<T>(pMin.x, pMax.y, pMin.z);
        else if(corner == 3) return Point3<T>(pMin.x, pMin.y, pMax.z);
        else if(corner == 4) return Point3<T>(pMin.x, pMax.y, pMax.z);
        else if(corner == 5) return Point3<T>(pMax.x, pMin.y, pMax.z);
        else if(corner == 6) return Point3<T>(pMax.x, pMax.y, pMin.z);
        else return Point3<T>(pMax.x, pMax.y, pMax.z);
    }
    Point3f Lerp(Point3f t) cosnt {
        return Point3f(Lerp(t.x, pMin.x, pMax.x),
                       Lerp(t.y, pMin.y, pMax.y),
                       Lerp(t.z, pMin.z, pMax.z));
    }
    Vector3f Offset(Point3f p) const {
        Vector3f o = p - pMin;
        if(pMax.x > pMin.x)
            o.x /= pMax.x - pMin.x;
        if(pMax.y > pMin.y)
            o.y /= pMax.y - pMin.y;
        if(pMax.z > pMin.z)
            o.z /= pMax.z - pMin.z;
        return o;
    }

    //几何性质辅助函数
    Vector3<T> Diagonal() const {
        return pMax - pMin;
    }

    T surfaceArea() const {
        Vector3<T> d = Diagonal();
        return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
    }

    T Volume() const {
        Vector3<T> d = Diagonal();
        return d.x * d.y * d.z;
    }

    int MaxLengthAxis() const {
        Vector3<T> d = Diagonal();
        if (d.x > d.y && d.x > d.z)
            return 0;
        else if (d.y > d.z)
            return 1;
        else
            return 2;
    }

    void BoundingSphere(Point3<T> *center, T *radius) const {
        *center = (pMin + pMax) / 2;
        *radius = Inside(*center, *this) ? Distance(*center, pMax) : 0;
    }

    bool isEmpty() const {
        return pMin.x >= pMax.x || pMin.y >= pMax.y || pMin.z >= pMax.z;
    }

    bool isDegenerate() const {
        return pMin.x > pMax.x || pMin.y > pMax.y || pMin.z > pMax.z;
    }

    template<typename U>
    explicit Bounds3(const Bounds3<U>& b) {
        if(b.isEmpty()){
            *this = Bounds3<T>();
        }
        else{
            pMin = Point3<T>(b.pMin);
            pMax = Point3<T>(b.pMax);
        }
    }

    bool operator==(const Bounds3<T>& b) cosnt {
        return b.pMin == pMin && b.pMax == pMax;
    }
    bool operator!=(const Bounds3<T> &b) const {
        return b.pMin != pMin || b.pMax != pMax;
    }
    //射线相交测试
    bool IntersectP(Point3f o, Vector3f d, Float tMax = Infinity, Float *hitt0 = nullptr,
                    Float *hitt1 = nullptr) const;
                    
    bool IntersectP(Point3f o, Vector3f d, Float tMax, Vector3f invDir,
                    const int dirIsNeg[3]) const;
};

template<typename T>
inline Bounds3<T> Union(const Bounds3<T>& b1, const Bounds3<T>& b2) {
    Bounds3<T> result;
    result.pMin = Min(b1.pMin, b2.pMin);
    result.pMax = Max(b1.pMax, b2.pMax);
    return result;
}
template <typename T>
inline Bounds3<T> Intersect(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    Bounds3<T> result;
    result.pMin = Max(b1.pMin, b2.pMin);
    result.pMax = Min(b1.pMax, b2.pMax);
    return result;
}
template <typename T>
inline bool Overlaps(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
    bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
    bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
    return (x && y && z);
}

template <typename T>
inline bool Inside(Point3<T> p, const Bounds3<T> &b) {
    return (p.x >= b.pMin.x && p.x <= b.pMax.x && p.y >= b.pMin.y && p.y <= b.pMax.y &&
            p.z >= b.pMin.z && p.z <= b.pMax.z);
}

template <typename T>
//Inside()的InsideExclusive()变体不认为位于上边界的点在包围盒内部
//它主要用于整数类型的包围盒
inline bool InsideExclusive(Point3<T> p, const Bounds3<T> &b) {
    return (p.x >= b.pMin.x && p.x < b.pMax.x && p.y >= b.pMin.y && p.y < b.pMax.y &&
            p.z >= b.pMin.z && p.z < b.pMax.z);
}
template <typename T, typename U>
inline auto DistanceSquared(Point3<T> p, const Bounds3<U> &b) {
    using TDist = decltype(T{} - U{});
    TDist dx = std::max<TDist>({0, b.pMin.x - p.x, p.x - b.pMax.x});
    TDist dy = std::max<TDist>({0, b.pMin.y - p.y, p.y - b.pMax.y});
    TDist dz = std::max<TDist>({0, b.pMin.z - p.z, p.z - b.pMax.z});
    return Sqr(dx) + Sqr(dy) + Sqr(dz);
}
template <typename T, typename U>
inline auto Distance(Point3<T> p, const Bounds3<U> &b) {
    auto dist2 = DistanceSquared(p, b);
    using TDist = typename TupleLength<decltype(dist2)>::type;
    return std::sqrt(TDist(dist2));
}

template <typename T, typename U>
inline Bounds3<T> Expand(const Bounds3<T> &b, U delta) {
    Bounds3<T> ret;
    ret.pMin = b.pMin - Vector3<T>(delta, delta, delta);
    ret.pMax = b.pMax + Vector3<T>(delta, delta, delta);
    return ret;
}



//记号
using Vector3f = Vector3<Float>;
using Point3f = Point3<Float>;
using Normal3f = Normal3<Float>;