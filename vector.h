#pragma once

#include "exceptions.h"
#include "matrix.h"
#include "myconcepts.h"

#include <cstddef>
#include <initializer_list>
#include <vector>

template<RingWithOne T>
class Vector {
public:
    Vector(size_t n) {
        data_.resize(n);
    }
    Vector(const Vector& other) = default;
    Vector(Vector&& other) = default;
    Vector& operator=(const Vector& other) = default;
    Vector& operator=(Vector&& other) = default;

    Vector(std::initializer_list<T> data) {
        data_.resize(data.size());
        size_t i = 0;
        for (auto it = data.begin(); it != data.end(); ++it, ++i) {
            data_[i] = *it;
        }
    }

    Vector(const Matrix<T>& mat) {
        if (mat.nsize() != 1 && mat.msize() != 1) {
            throw WrongSizeException();
        }
        if (mat.msize() == 1) {
            data_.resize(mat.nsize());
            for (size_t i = 0; i < data_.size(); ++i) {
                data_[i] = mat[i][0];
            }
        } else {
            data_.resize(mat.msize());
            for (size_t i = 0; i < data_.size(); ++i) {
                data_[i] = mat[0][i];
            }
        }
    }

    size_t size() const {
        return data_.size();
    }

    const T& operator[](size_t i) const {
        return data_[i];
    }
    T& operator[](size_t i) {
        return data_[i];
    }

    Vector operator+(const Vector& other) const {
        if (size() != other.size()) {
            throw WrongSizeException();
        }
        Vector ans;
        for (size_t i = 0; i < size(); ++i) {
            ans[i] = other[i] + data_[i];
        }
        return ans;
    }
    Vector operator-(const Vector& other) const {
        if (size() != other.size()) {
            throw WrongSizeException();
        }
        Vector ans(size());
        for (size_t i = 0; i < size(); ++i) {
            ans[i] = data_[i] - other[i];
        }
        return ans;
    }
    T operator*(const Vector& other) const {
        if (size() != other.size()) {
            throw WrongSizeException();
        }
        T ans = T::ZERO();
        for (size_t i = 0; i < size(); ++i) {
            ans += data_[i] * other[i];
        }
        return ans;
    }
    Vector operator*(const T& lambda) const {
        Vector ans = *this;
        for (size_t i = 0; i < size(); ++i) {
            ans[i] *= lambda;
        }
        return ans;
    }
    Vector operator-() const {
        Vector ans(size());
        for (size_t i = 0; i < size(); ++i) {
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

    operator Matrix<T>() {
        Matrix<T> ans(size(), 1);
        for (size_t i = 0; i < size(); ++i) {
            ans[i][0] = data_[i];
        }
        return ans;
    }

    friend Vector operator*(const T& lambda, const Vector& vec) {
        return vec * lambda;
    }

    friend std::ostream &operator<<(std::ostream& out, const Vector& vec) {
        out << "(";
        for (size_t i = 0; i < vec.size(); ++i) {
            out << (i ? ", " : "") << vec[i];
        }
        out << ")";
        return out;
    }

private:
    std::vector<T> data_;
};
