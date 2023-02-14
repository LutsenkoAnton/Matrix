#pragma once

#include "matrix.h"
#include "myconcepts.h"

#include <array>
#include <cstddef>
#include <initializer_list>

template<RingWithOne T, size_t n>
class Vector {
public:
    Vector() = default;
    Vector(const Vector& other) = default;
    Vector(Vector&& other) = default;
    Vector& operator=(const Vector& other) = default;
    Vector& operator=(Vector&& other) = default;

    Vector(std::initializer_list<T> data) {
        size_t i = 0;
        for (auto it = data.begin(); it != data.end(); ++it, ++i) {
            data_[i] = *it;
        }
    }

    const T& operator[](size_t i) const {
        return data_[i];
    }
    T& operator[](size_t i) {
        return data_[i];
    }

    Vector operator+(const Vector& other) const {
        Vector ans;
        for (size_t i = 0; i < n; ++i) {
            ans[i] = other[i] + data_[i];
        }
        return ans;
    }
    Vector operator-(const Vector& other) const {
        Vector ans;
        for (size_t i = 0; i < n; ++i) {
            ans[i] = other[i] - data_[i];
        }
        return ans;
    }
    T operator*(const Vector& other) const {
        T ans = T::ZERO;
        for (size_t i = 0; i < n; ++i) {
            ans += data_[i] * other[i];
        }
        return ans;
    }
    Vector operator*(const T& lambda) const {
        Vector ans = *this;
        for (size_t i = 0; i < n; ++i) {
            ans[i] *= lambda;
        }
        return ans;
    }
    Vector operator-() const {
        Vector ans;
        for (size_t i = 0; i < n; ++i) {
            ans[i] = -data_[i];
        }
        return ans;
    }

    Vector& operator+=(const Vector& other) {
        return *this = *this + other;
    }
    Vector& operator-=(const Vector& other) {
        return *this = *this - other;
    }
    Vector& operator*=(const T& lambda) {
        return *this = *this * lambda;
    }

    operator Matrix<T, n, 1>() {
        Matrix<T, n, 1> ans;
        for (size_t i = 0; i < n; ++i) {
            ans[i][0] = data_[i];
        }
        return ans;
    }

    friend Vector operator*(const T& lambda, const Vector& vec) {
        return vec * lambda;
    }

    friend std::ostream &operator<<(std::ostream& out, const Vector& vec) {
        out << "(";
        for (size_t i = 0; i < n; ++i) {
            out << (i ? ", " : "") << vec[i];
        }
        out << ")";
        return out;
    }

private:
    std::array<T, n> data_;
};
