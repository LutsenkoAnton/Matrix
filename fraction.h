#pragma once

#include <compare>
#include <iostream>
#include <numeric>

class Fraction {
public:
    Fraction();
    Fraction(int a);  // NOLINT[google-explicit-constructor]
    Fraction(int a, int b);

    bool operator==(const Fraction& other) const;
    bool operator!=(const Fraction& other) const;
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

    int numerator_;
    int denominator_;
};
