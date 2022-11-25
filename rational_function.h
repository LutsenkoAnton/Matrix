#pragma once

#include "poly.h"

#include <exception>
#include <iostream>

class NotDivisibleException: public std::exception {
public:
    const char* what() const noexcept override {
        return "The numerator is not divisible by the denominator";
    }
};

template<typename T>
class RationalFunction {
public:
    RationalFunction(): numerator(T(0)), denominator(T(1)) {}
    RationalFunction(const RationalFunction& other) = default;
    RationalFunction(RationalFunction&& other) = default;
    RationalFunction& operator=(const RationalFunction& other) = default;
    RationalFunction& operator=(RationalFunction&& other) = default;

    RationalFunction(const T& a): numerator(a), denominator(T(1)) {}
    RationalFunction(const Poly<T>& a): numerator(a), denominator(T(1)) {}
    RationalFunction(const Poly<T>& a, const Poly<T>& b): numerator(a), denominator(b) {}

    bool operator==(const RationalFunction& other) const {
        return numerator * other.denominator == denominator * other.numerator;
    }
    bool operator!=(const RationalFunction& other) const {
        return !operator==(other);
    }

    RationalFunction operator+(const RationalFunction& other) const {
        auto ans = RationalFunction(numerator * other.denominator + denominator * other.numerator, denominator * other.denominator);
        ans.try_devide();
        return ans;
    }
    RationalFunction operator-(const RationalFunction& other) const {
        return *this + (-other);
    }
    RationalFunction operator*(const RationalFunction& other) const {
        auto ans = RationalFunction(numerator * other.numerator, denominator * other.denominator);
        ans.try_devide();
        return ans;
    }
    RationalFunction operator/(const RationalFunction& other) const {
        auto ans = RationalFunction(numerator * other.denominator, denominator * other.numerator);
        ans.try_devide();
        return ans;
    }
    RationalFunction operator-() const {
        return RationalFunction(-numerator, denominator);
    }

    RationalFunction& operator+=(const RationalFunction& other) {
        return *this = *this + other;
    }
    RationalFunction& operator-=(const RationalFunction& other) {
        return *this = *this - other;
    }
    RationalFunction& operator*=(const RationalFunction& other) {
        return *this = *this * other;
    }
    RationalFunction& operator/=(const RationalFunction& other) {
        return *this = *this / other;
    }

    void try_devide() {
        if (numerator == Poly<T>(0)) {
            denominator = Poly<T>(1);
            return;
        }
        size_t n = numerator.deg();
        size_t m = denominator.deg();
        Poly<T> curans;
        T am = denominator.SeniorCoefficient();
        auto cpnumerator = numerator;
        for (int i = n - m; i >= 0; --i) {
            T an = cpnumerator.SeniorCoefficient();
            curans += Poly<T>(an / am, i);
            cpnumerator -= denominator * Poly<T>(an / am, i);
        }
        cpnumerator.Clean();
        if (cpnumerator == Poly<T>(T(0))) {
            numerator = curans;
            denominator = Poly<T>(T(1));
        }
    }

    Poly<T> force_devide() {
        try_devide();
        if (denominator != Poly<T>(T(1))) {
            throw NotDivisibleException();
        }
        return numerator;
    }

    friend std::ostream& operator<<(std::ostream& stream, const RationalFunction p) {
        stream << '(' << p.numerator << ") / (" << p.denominator << ")";
        return stream;
    }

private:
    Poly<T> numerator;
    Poly<T> denominator;
};


