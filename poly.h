#pragma once

#include "mymath.h"

#include <functional>
#include <initializer_list>
#include <iostream>
#include <map>
#include <vector>

template <typename T>
class Poly {
public:
    Poly() {}
    Poly(const Poly& other) = default;
    Poly(const Poly&& other) = default;
    Poly& operator=(const Poly& other) = default;
    Poly& operator=(const Poly&& other) = default;
    Poly(const T& coefficient): Poly(coefficient, 0) {
        Clean();
    }
    Poly(const T& coefficient, size_t power): coefficients_({{power, coefficient}}) {}
    Poly(const std::vector<T>& coefficients) {
        for (size_t power = 0; power < coefficients.size(); ++power) {
            if (coefficients[power] != 0) {
                coefficients_[power] = coefficients[power];
            }
        }
    }
    Poly(const std::map<size_t, T, std::greater<size_t>> coefficients): coefficients_(coefficients) {
        Clean();
    }
    Poly(const std::initializer_list<T>& coefficients): Poly(std::vector<T>(coefficients)) {}

    T operator() (const T& x) const {
        T sum = 0;
        for (const auto& [power, coefficient] : coefficients_) {
            sum += fastpow(x, power) * coefficient;
        }
        return sum;
    }

    bool operator==(const Poly& other) const = default;
    bool operator!=(const Poly& other) const = default;

    Poly operator+(const Poly& other) const {
        Poly ans(*this);
        for (const auto& [power, coefficient] : other.coefficients_) {
            ans.coefficients_[power] += coefficient;
        }
        ans.Clean();
        return ans;
    }
    Poly operator-() const {
        Poly ans(*this);
        for (auto &[power, coefficient] : ans.coefficients_) {
            coefficient = -coefficient;
        }
        return ans;
    }
    Poly operator-(const Poly& other) const {
        return *this + -other;
    }
    Poly operator*(const Poly& other) const {
        Poly ans;
        for (const auto& [power1, coefficient1] : coefficients_) {
            for (const auto& [power2, coefficient2] : other.coefficients_) {
                ans.coefficients_[power1 + power2] += coefficient1 * coefficient2;
            }
        }
        ans.Clean();
        return ans;
    }

    Poly& operator+=(const Poly& other) {
        *this = *this + other;
        return *this;
    }
    Poly& operator-=(const Poly& other) {
        *this = *this - other;
        return *this;
    }
    Poly& operator*=(const Poly& other) {
        *this = *this * other;
        return *this;
    }

    std::map<size_t, T, std::greater<size_t>> GetCoefficients() const {
        return coefficients_;
    }

    std::map<size_t, T, std::greater<size_t>>& GetCoefficients() {
        return coefficients_;
    }

    void Clean() {
        for (auto it = coefficients_.begin(); it != coefficients_.end();) {
            if (it -> second == T(0)) {
                it = coefficients_.erase(it);
            } else {
                ++it;
            }
        }
    }

    T SeniorCoefficient() const {
        if (coefficients_.empty()) return 0;
        return coefficients_.begin() -> second;
    }

    size_t deg() const {
        if (coefficients_.empty()) return 0;
        return coefficients_.begin() -> first;
    }

    friend std::ostream& operator<<(std::ostream& stream, const Poly& polynomial) {
        if (polynomial.coefficients_.empty()) {
            stream << "0";
            return stream;
        }
        bool first = true;
        for (const auto& [power, coefficient] : polynomial.coefficients_) {
            if (coefficient == T(0)) continue;
            if (first) {
                first = false;
                if (coefficient != T(1) || power == 0) {
                    stream << coefficient;
                }
                if (power == 0) {
                    continue;
                }
                if (power == 1) {
                    stream << "x";
                    continue;
                }
                stream << "x^" << power;
                continue;
            }
            stream << ' ' << (coefficient < T(0) ? '-' : '+') << ' ';
            if (myabs(coefficient) != T(1) || power == 0) {
                stream << myabs(coefficient);
            }
            if (power == 0) {
                continue;
            }
            stream << 'x';
            if (power == 1) {
                continue;
            }
            stream << '^' << power;
        }
        return stream;
    }

private:
    std::map<size_t, T, std::greater<size_t>> coefficients_;
};

