#pragma once

#include <compare>
#include <iostream>

#include "mymath.h"

class Float {
public:
    Float() = default;
    Float(const Float& other) = default;
    Float(Float&& other) = default;
    Float& operator=(const Float& other) {
        data_ = other.data_;
        return *this;
    }
    Float& operator=(Float&& other) {
        data_ = other.data_;
        return *this;
    }

    Float(double num): data_(num) {
        if (abs(num) < EPS) { // Gets rid of -0s
            data_ = 0;
        }
    }

    bool operator==(const Float& other) const {
        return abs(other.data_ - data_) < EPS;
    }
    bool operator!=(const Float& other) const {
        return !operator==(other);
    }
    std::strong_ordering operator<=>(const Float& other) const {
        if (*this == other) {
            return std::strong_ordering::equal;
        }
        if (data_ < other.data_) {
            return std::strong_ordering::less;
        }
        return std::strong_ordering::greater;
    }

    Float operator+(const Float& other) const {
        return Float(data_ + other.data_);
    }
    Float operator-(const Float& other) const {
        return Float(data_ - other.data_);
    }
    Float operator*(const Float& other) const {
        return Float(data_ * other.data_);
    }
    Float operator/(const Float& other) const {
        return Float(data_ / other.data_);
    }
    Float operator-() const {
        return Float(-data_);
    }
    Float& operator+=(const Float& other) {
        return *this = *this + other;
    }
    Float& operator-=(const Float& other) {
        return *this = *this - other;
    }
    Float& operator*=(const Float& other) {
        return *this = *this * other;
    }
    Float& operator/=(const Float& other) {
        return *this = *this / other;
    }

    friend std::ostream& operator<<(std::ostream& out, const Float& i) {
        out << i.data_;
        return out;
    }

    friend std::istream& operator<<(std::istream& in, Float& i) {
        in >> i.data_;
        return in;
    }

    static Float ONE() {
        return Float(1);
    }
    static Float ZERO() {
        return Float(0);
    }

private:
    const double EPS = 1e-6;
    double data_;
};

inline Float operator "" _f(long double i) {
    return Float(static_cast<double>(i));
}

inline Float operator "" _f(unsigned long long i) {
    return Float(static_cast<double>(i));
}
