#pragma once

#include <cmath>
#include <compare>
#include <iostream>
#include <numeric>

class Fraction {
public:
    Fraction();
    Fraction(int64_t a);
    Fraction(int64_t a, int64_t b);
    Fraction(const Fraction& other) = default;
    Fraction(Fraction&& other) = default;
    Fraction& operator=(const Fraction& other) = default;
    Fraction& operator=(Fraction&& other) = default;

    bool operator==(const Fraction& other) const = default;
    bool operator!=(const Fraction& other) const = default;
    std::strong_ordering operator<=>(const Fraction& other) const;

    Fraction operator+(const Fraction& other) const;
    Fraction operator-(const Fraction& other) const;
    Fraction operator*(const Fraction& other) const;
    Fraction operator/(const Fraction& other) const;
    Fraction operator-() const;

    Fraction& operator+=(const Fraction& other);
    Fraction& operator-=(const Fraction& other);
    Fraction& operator*=(const Fraction& other);
    Fraction& operator/=(const Fraction& other);

    friend std::ostream& operator<<(std::ostream& stream, const Fraction& fraction);

private:
    void Shorten();

    int64_t numerator_;
    int64_t denominator_;
};
