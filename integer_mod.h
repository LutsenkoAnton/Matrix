#pragma once

#include <compare>
#include <iostream>

#include "mymath.h"

class IntegerMod {
public:
    IntegerMod() = default;
    IntegerMod(const IntegerMod& other) = default;
    IntegerMod(IntegerMod&& other) = default;
    IntegerMod& operator=(const IntegerMod& other) = default;
    IntegerMod& operator=(IntegerMod&& other) = default;

    IntegerMod(int num): data_(num) {
        Normalize();
    }

    IntegerMod(long long num): data_(num) {
        Normalize();
    }

    bool operator==(const IntegerMod& other) const = default;
    bool operator!=(const IntegerMod& other) const = default;
    std::strong_ordering operator<=>(const IntegerMod& other) const = default;

    IntegerMod operator+(const IntegerMod& other) const {
        return IntegerMod(data_ + other.data_);
    }
    IntegerMod operator-(const IntegerMod& other) const {
        return IntegerMod(data_ - other.data_);
    }
    IntegerMod operator*(const IntegerMod& other) const {
        return IntegerMod(data_ * other.data_);
    }
    IntegerMod operator/(const IntegerMod& other) const {
        return *this * fastpow(other, MOD - 2);
    }
    IntegerMod operator-() const {
        return IntegerMod(-data_);
    }
    IntegerMod& operator+=(const IntegerMod& other) {
        return *this = *this + other;
    }
    IntegerMod& operator-=(const IntegerMod& other) {
        return *this = *this - other;
    }
    IntegerMod& operator*=(const IntegerMod& other) {
        return *this = *this * other;
    }
    IntegerMod& operator/=(const IntegerMod& other) {
        return *this = *this / other;
    }

    friend std::ostream& operator<<(std::ostream& out, const IntegerMod& i) {
        out << i.data_;
        return out;
    }

    friend std::istream& operator<<(std::istream& in, IntegerMod& i) {
        in >> i.data_;
        return in;
    }

    static IntegerMod ONE() {
        return IntegerMod(1);
    }
    static IntegerMod ZERO() {
        return IntegerMod(0);
    }

private:
    static const long long MOD = 1e9 + 9;
    void Normalize() {
        data_ %= MOD;
        while (data_ > MOD / 2) {
            data_ -= MOD;
        }
        while (data_ < -(MOD / 2)) {
            data_ += MOD;
        }
    }
    long long data_;
};

inline IntegerMod operator "" _im(unsigned long long i) {
    return IntegerMod(static_cast<long long>(i));
}
