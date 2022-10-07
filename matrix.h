#pragma once

#include <array>
#include <initializer_list>
#include <iostream>
#include <numeric>

#include "fraction.h"

template <size_t n, size_t m>
class Matrix {
public:
    Matrix();
    Matrix(const Matrix& other);
    explicit Matrix(std::array<std::array<Fraction, m>, n> data);
    Matrix(std::initializer_list<std::initializer_list<Fraction>> data);
    Matrix& operator=(const Matrix& other);

    Matrix operator+(const Matrix& other) const;

    Matrix& operator+=(const Matrix& other);

    template <size_t k>
    Matrix<n, k> operator*(const Matrix<m, k>& other) const;

    Matrix operator*(const Fraction lambda) const;

    Matrix& operator*=(const Fraction lambda);

    template <size_t k>
    Matrix<n, m + k> operator|(const Matrix<n, k>& other) const;

    std::array<Fraction, m>& operator[](size_t pos);
    std::array<Fraction, m> operator[](size_t pos) const;

    void Gauss();

    Matrix<m, n> Transpose() const;

    template <size_t x, size_t y>
    void Slice(Matrix<x, y>& to, size_t start_i = 0, size_t start_j = 0) const;

    Fraction Trace() const;

private:
    std::array<std::array<Fraction, m>, n> data_;
};

template <size_t n, size_t m>
Matrix<n, m> operator*(const Fraction lambda, const Matrix<n, m>& matrix);

template <size_t n, size_t m>
std::ostream& operator<<(std::ostream& stream, const Matrix<n, m>& matrix);

template <size_t n>
class SquareMatrix : public Matrix<n, n> {
public:
    SquareMatrix();
    SquareMatrix(Fraction lambda);  // NOLINT[google-explicit-constructor]

    SquareMatrix Power(size_t indicator) const;

    SquareMatrix Inverse() const;

    static SquareMatrix IdentityMatrix();
};
