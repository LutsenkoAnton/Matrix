#pragma once

#include <initializer_list>
#include <iostream>
#include <vector>

template <typename T>
class Poly {
public:
    Poly() {}
    Poly(const T& coefficient): coefficients_({coefficient}) {}
    Poly(const std::vector<T>& coefficients): coefficients_(coefficients) {}
    Poly(const std::initializer_list<T>& coefficients): coefficients_(coefficients) {}
    Poly& operator=(const Poly& other) = default;

    T operator() (const T& x) const {
        T power = 1;
        T sum = 0;
        for (size_t i = 0; i < coefficients_.size(); ++i) {
            sum += power * coefficients_[i];
            power *= x;
        }
        return sum;
    }

    bool operator==(const Poly& other) const = default;
    bool operator!=(const Poly& other) const = default;

    Poly operator+(const Poly& other) const {
        Poly ans(std::vector<T>(std::max(coefficients_.size(), other.coefficients_.size())));
        for (size_t i = 0; i < coefficients_.size(); ++i) {
            ans.coefficients_[i] = coefficients_[i];
        }
        for (size_t i = 0; i < other.coefficients_.size(); ++i) {
            ans.coefficients_[i] += other.coefficients_[i];
        }
        return ans;
    }
    Poly operator-() const {
        Poly ans(coefficients_);
        for (size_t i = 0; i < coefficients_.size(); ++i) {
            ans.coefficients_[i] = -ans.coefficients_[i];
        }
        return ans;
    }
    Poly operator-(const Poly& other) const {
        return *this + -other;
    }
    Poly operator*(const Poly& other) const {
        Poly ans(std::vector<T>(coefficients_.size() + other.coefficients_.size() - 1));
        for (size_t i = 0; i < coefficients_.size(); ++i) {
            for (size_t j = 0; j < other.coefficients_.size(); ++j) {
                ans.coefficients_[i + j] += coefficients_[i] * other.coefficients_[j];
            }
        }
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

    friend std::ostream& operator<<(std::ostream& stream, const Poly& polynomial) {
        if (polynomial.coefficients_.empty()) {
            stream << "0";
            return stream;
        }
        bool first = true;
        for (size_t p = 0; p < polynomial.coefficients_.size(); ++p) {
            size_t power = polynomial.coefficients_.size() - 1 - p;
            auto coefficient = polynomial.coefficients_[power];
            if (coefficient == 0) continue;
            if (first) {
                first = false;
                if (coefficient != 1) {
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
            stream << ' ' << (coefficient < 0 ? '-' : '+') << ' ' << std::abs(coefficient);
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
        std::vector<T> coefficients_;
};

