#pragma once

#include "integer.h"
#include "myconcepts.h"

#include <compare>
#include <iostream>

template<EuclideanRing T>
class Fraction {
public:
    Fraction() : numerator_(T::ZERO), denominator_(T::ONE) {}
    Fraction(T a) : numerator_(a), denominator_(T::ONE) {}
    Fraction(T a, T b) : numerator_(a), denominator_(b) {
        Shorten();
    }
    Fraction(const Fraction& other) = default;
    Fraction(Fraction&& other) = default;
    Fraction& operator=(const Fraction& other) = default;
    Fraction& operator=(Fraction&& other) = default;

    bool operator==(const Fraction& other) const {
        return other.denominator_ * numerator_ == other.numerator_ * denominator_;
    }
    bool operator!=(const Fraction& other) const {
        return !operator==(other);
    }
    std::strong_ordering operator<=>(const Fraction& other) const {
        return other.denominator_ * numerator_ <=> other.numerator_ * denominator_;
    }

    Fraction operator+(const Fraction& other) const {
        return Fraction(numerator_ * other.denominator_ + denominator_ * other.numerator_, denominator_ * other.denominator_);
    }
    Fraction operator-(const Fraction& other) const {
        return Fraction(numerator_ * other.denominator_ - denominator_ * other.numerator_, denominator_ * other.denominator_);
    }
    Fraction operator*(const Fraction& other) const {
        return Fraction(numerator_ * other.numerator_, denominator_ * other.denominator_);
    }
    Fraction operator/(const Fraction& other) const {
        return Fraction(numerator_ * other.denominator_, denominator_ * other.numerator_);
    }
    Fraction operator-() const {
        return Fraction(-numerator_, denominator_);
    }

    Fraction& operator+=(const Fraction& other) {
        return *this = *this + other;
    }
    Fraction& operator-=(const Fraction& other) {
        return *this = *this - other;
    }
    Fraction& operator*=(const Fraction& other) {
        return *this = *this * other;
    }
    Fraction& operator/=(const Fraction& other) {
        return *this = *this / other;
    }

    T GetIntegerPart() const {
        return numerator_ / denominator_;
    }
    Fraction GetFractionalPart() const {
        return Fraction(numerator_ % denominator_, denominator_);
    }

    friend std::ostream& operator<<(std::ostream& stream, const Fraction& fraction) {
        if (fraction.denominator_ == T::ONE) {
            stream << fraction.numerator_;
            return stream;
        }
        stream << "(" << fraction.numerator_ << ")/(" << fraction.denominator_ << ")";
        return stream;
    }

    static inline const Fraction ZERO = Fraction();
    static inline const Fraction ONE = Fraction(T::ONE);

private:
    void Shorten() {
        T d = gcd(numerator_, denominator_);
        numerator_ /= d;
        denominator_ /= d;
    }

    T numerator_;
    T denominator_;
};

Fraction<Integer> operator "" _fi(unsigned long long a) {
    return Fraction<Integer>(Integer(static_cast<long long>(a)));
}
