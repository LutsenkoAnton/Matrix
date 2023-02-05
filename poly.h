#pragma once

#include "mymath.h"
#include "myconcepts.h"

#include <functional>
#include <initializer_list>
#include <iostream>
#include <map>
#include <vector>

template <RingWithOne T>
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
            if (coefficients[power] != T::ZERO) {
                coefficients_[power] = coefficients[power];
            }
        }
    }
    Poly(const std::map<size_t, T, std::greater<size_t>> coefficients): coefficients_(coefficients) {
        Clean();
    }
    Poly(const std::initializer_list<T>& coefficients): Poly(std::vector<T>(coefficients)) {}

    T operator() (const T& x) const {
        T sum = T::ZERO;
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

    Poly operator/(const Poly& other) const requires Field<T> {
        if (*this == Poly::ZERO) {
            return Poly::ZERO;
        }
        size_t n = deg();
        size_t m = other.deg();
        Poly<T> curans;
        T am = other.SeniorCoefficient();
        auto cpthis = *this;
        for (int i = n - m; i >= 0; --i) {
            T an = cpthis.SeniorCoefficient();
            curans += Poly(an / am, i);
            cpthis -= other * Poly(an / am, i);
        }
        curans.Clean();
        return curans;
    }

    Poly operator%(const Poly& other) const requires Field<T> {
        if (*this == Poly::ZERO) {
            return Poly::ZERO;
        }
        size_t n = deg();
        size_t m = other.deg();
        Poly<T> curans;
        T am = other.SeniorCoefficient();
        auto cpthis = *this;
        for (int i = n - m; i >= 0; --i) {
            T an = cpthis.SeniorCoefficient();
            curans += Poly(an / am, i);
            cpthis -= other * Poly(an / am, i);
        }
        cpthis.Clean();
        return cpthis;
    }

    Poly& operator+=(const Poly& other) {
        return *this = *this + other;
    }
    Poly& operator-=(const Poly& other) {
        return *this = *this - other;
    }
    Poly& operator*=(const Poly& other) {
        return *this = *this * other;
    }
    Poly& operator/=(const Poly& other) requires Field<T> {
        return *this = *this / other;
    }
    Poly& operator%=(const Poly& other) requires Field<T> {
        return *this = *this % other;
    }

    std::map<size_t, T, std::greater<size_t>> GetCoefficients() const {
        return coefficients_;
    }

    std::map<size_t, T, std::greater<size_t>>& GetCoefficients() {
        return coefficients_;
    }

    void Clean() {
        for (auto it = coefficients_.begin(); it != coefficients_.end();) {
            if (it -> second == T::ZERO) {
                it = coefficients_.erase(it);
            } else {
                ++it;
            }
        }
    }

    T SeniorCoefficient() const {
        if (coefficients_.empty()) return T::ZERO;
        return coefficients_.begin() -> second;
    }

    size_t deg() const {
        if (coefficients_.empty()) return 0;
        return coefficients_.begin() -> first;
    }

    friend std::ostream& operator<<(std::ostream& stream, const Poly& polynomial) {
        if (polynomial.coefficients_.empty()) {
            stream << T::ZERO;
            return stream;
        }
        bool first = true;
        for (const auto& [power, coefficient] : polynomial.coefficients_) {
            if (coefficient == T::ZERO) continue;
            if (first) {
                first = false;
                if (coefficient != T::ONE || power == 0) {
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
            stream << ' ' << (coefficient < T::ZERO ? '-' : '+') << ' ';
            if (abs(coefficient) != T::ONE || power == 0) {
                stream << abs(coefficient);
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

    static inline const Poly ZERO = Poly(T::ZERO);
    static inline const Poly ONE = Poly(T::ONE);

private:
    std::map<size_t, T, std::greater<size_t>> coefficients_;
};

