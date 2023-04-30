#pragma once

#include <compare>
#include <iostream>

#include "mymath.h"

class Integer {
public:
    Integer() = default;
    Integer(const Integer& other) = default;
    Integer(Integer&& other) = default;
    Integer& operator=(const Integer& other) = default;
    Integer& operator=(Integer&& other) = default;

    Integer(int num): data_(num) {}
    Integer(long long num): data_(num) {}

    bool operator==(const Integer& other) const = default;
    bool operator!=(const Integer& other) const = default;
    std::strong_ordering operator<=>(const Integer& other) const = default;

    Integer operator+(const Integer& other) const {
        return Integer(data_ + other.data_);
    }
    Integer operator-(const Integer& other) const {
        return Integer(data_ - other.data_);
    }
    Integer operator*(const Integer& other) const {
        return Integer(data_ * other.data_);
    }
    Integer operator/(const Integer& other) const {
        return Integer(data_ / other.data_);
    }
    Integer operator%(const Integer& other) const {
        return Integer(data_ % other.data_);
    }
    Integer operator-() const {
        return Integer(-data_);
    }
    Integer& operator+=(const Integer& other) {
        return *this = *this + other;
    }
    Integer& operator-=(const Integer& other) {
        return *this = *this - other;
    }
    Integer& operator*=(const Integer& other) {
        return *this = *this * other;
    }
    Integer& operator/=(const Integer& other) {
        return *this = *this / other;
    }
    Integer& operator%=(const Integer& other) {
        return *this = *this % other;
    }

    friend std::ostream& operator<<(std::ostream& out, const Integer& i) {
        out << i.data_;
        return out;
    }

    friend std::istream& operator<<(std::istream& in, Integer& i) {
        in >> i.data_;
        return in;
    }

    static Integer ONE() {
        return Integer(1);
    }
    static Integer ZERO() {
        return Integer(0);
    }

private:
    long long data_;
};

inline Integer operator "" _i(unsigned long long i) {
    return Integer(static_cast<long long>(i));
}
