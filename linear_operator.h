#pragma once

#include "matrix.h"
#include "myconcepts.h"
#include "poly.h"
#include "vector.h"
#include "vector_space.h"

#include <cstddef>
#include <iterator>

template<Field T, size_t n>
class LinearOperator {
public:
    LinearOperator(const Matrix<T, n>& data): data_(data) {}
    LinearOperator(const VectorSpace<T, n>& im, const VectorSpace<T, n>& ker) : LinearOperator(im.GetBasisMatrix() * ker.GetEquasion()) {}
    LinearOperator operator+(const LinearOperator& other) const {
        return LinearOperator(data_ + other.data_);
    }
    LinearOperator operator-(const LinearOperator& other) const {
        return LinearOperator(data_ - other.data_);
    }
    LinearOperator operator*(const LinearOperator& other) const {
        return LinearOperator(data_ * other.data_);
    }
    LinearOperator operator*(const T& lambda) const {
        return LinearOperator(data_ * lambda);
    }
    friend LinearOperator operator*(const T& lambda, const LinearOperator& op) {
        return LinearOperator(op.data_ * lambda);
    }
    LinearOperator& operator+=(const LinearOperator& other) {
        return *this = *this + other;
    }
    LinearOperator& operator-=(const LinearOperator& other) {
        return *this = *this - other;
    }
    LinearOperator& operator*=(const LinearOperator& other) {
        return *this = *this * other;
    }
    LinearOperator& operator*=(const T& lambda) {
        return *this = *this * lambda;
    }
    LinearOperator Power(size_t indicator) const {
        return LinearOperator(data_.Power(indicator));
    }

    friend T Det(const LinearOperator& op) {
        return Det(op.data_);
    }

    Poly<T> CharPoly() const {
        return data_.CharPoly();
    }

    T Trace() const {
        return data_.Trace();
    }

    VectorSpace<T, n> Im() const {
        std::array<Vector<T, n>, n> vectors;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                vectors[i][j] = data_[j][i];
            }
        }
        return VectorSpace<T, n>(vectors);
    }

    VectorSpace<T, n> ker() const {
        return VectorSpace<T, n>(data_);
    }

    Matrix<T, n> GetData() const {
        return data_;
    }

    static inline const LinearOperator ZERO = LinearOperator();
    static inline const LinearOperator ONE = LinearOperator(Matrix<T, n>::IdentityMatrix());
private:
    Matrix<T, n> data_;
};
