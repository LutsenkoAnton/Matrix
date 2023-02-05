#pragma once

#include "fraction.h"
#include "permutation.h"
#include "poly.h"
#include "myconcepts.h"

#include <array>
#include <cassert>
#include <initializer_list>
#include <iostream>
#include <exception>
#include <numeric>
#include <type_traits>

class WrongSizeException: public std::exception {
public:
    const char* what() const noexcept override {
        return "Wrong sizes of matrixes";
    }
};

class SingularMatrixException: public std::exception {
public:
    const char* what() const noexcept override {
        return "Matrix is singular";
    }
};

template <RingWithOne T, size_t n, size_t m = n>
class Matrix {
public:
    Matrix() = default;
    Matrix(const Matrix& other) = default;
    Matrix(Matrix&& other) = default;
    Matrix& operator=(const Matrix& other) = default;
    Matrix& operator=(Matrix&& other) = default;

    Matrix(std::initializer_list<std::initializer_list<T>> data) {
        size_t i = 0;
        for (const auto& row : data) {
            size_t j = 0;
            for (const auto& elem : row) {
                data_[i][j] = elem;
                ++j;
            }
            ++i;
        }
        for (size_t i = 0; i < data_.size(); ++i) {
            if (data_[i].size() != data_[0].size()) throw WrongSizeException();
        }
    }

    bool operator==(const Matrix& other) const = default;
    bool operator!=(const Matrix& other) const = default;

    std::pair<size_t, size_t> size() const {
        return std::make_pair(n, m);
    }

    Matrix operator+(const Matrix& other) const {
        Matrix ans;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                ans[i][j] = data_[i][j] + other[i][j];
            }
        }
        return ans;
    }

    Matrix& operator+=(const Matrix& other) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                data_[i][j] += other[i][j];
            }
        }
        return *this;
    }

    Matrix operator-() const {
        Matrix ans = *this;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                ans[i][j] = -ans[i][j];
            }
        }
        return ans;
    }

    Matrix operator-(const Matrix& other) const {
        return *this + -other;
    }

    Matrix& operator-=(const Matrix& other) {
        return *this = *this - other;
    }

    template<size_t k>
    Matrix<T, n, k> operator*(const Matrix<T, m, k>& other) const {
        Matrix<T, n, k> a;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                for (size_t l = 0; l < k; ++l) {
                    a[i][l] += data_[i][j] * other[j][l];
                }
            }
        }
        return a;
    }

    Matrix operator*(const T& lambda) const {
        Matrix ans;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                ans[i][j] = lambda * data_[i][j];
            }
        }
        return ans;
    }

    Matrix& operator*=(const T& lambda) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                data_[i][j] *= lambda;
            }
        }
    }

    // writes other to the right from *this
    template<size_t k>
    Matrix<T, n, m + k> operator|(const Matrix<T, n, k>& other) const {
        Matrix<T, n, m + k> a;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                a[i][j] = data_[i][j];
            }
            for (size_t j = 0; j < k; ++j) {
                a[i][j + m] = other[i][j];
            }
        }
        return a;
    }

    std::array<T, m>& operator[](size_t pos) {
        return data_[pos];
    }
    const std::array<T, m>& operator[](size_t pos) const {
        return data_[pos];
    }

    void Gauss() {
        size_t start = 0;
        for (size_t j = 0; j < m; ++j) {
            size_t with_non_zero_coefficient = 0;
            bool found = false;
            for (size_t i = start; i < n; ++i) {
                if (data_[i][j] != T::ZERO) {
                    found = true;
                with_non_zero_coefficient = i;
                    break;
                }
            }
            if (!found) {
                continue;
            }
            swap(data_[start], data_[with_non_zero_coefficient]);
            for (size_t i = j + 1; i < m; ++i) {
                data_[start][i] /= data_[start][j];
            }
            data_[start][j] = T::ONE;
            for (size_t i = 0; i < n; ++i) {
                if (i == start) {
                    continue;
                }
                if (data_[i][j] == T::ZERO) {
                    continue;
                }
                auto lambda = -data_[i][j];
                for (size_t k = j; k < m; ++k) {
                    data_[i][k] += lambda * data_[start][k];
                }
            }
            ++start;
        }
    }

    Matrix Transpose() const {
        Matrix<T, m, n> ans;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                ans[j][i] = data_[i][j];
            }
        }
        return ans;
    }

    template<size_t len_i, size_t len_j>
    Matrix<T, len_i, len_j> Slice(size_t start_i = 0, size_t start_j = 0) const {
        if (len_i + start_i > n || len_j + start_j > m) throw WrongSizeException();
        Matrix<T, len_i, len_j> ans;
        for (size_t i = 0; i < len_i; ++i) {
            for (size_t j = 0; j < len_j; ++j) {
                ans[i][j] = data_[i + start_i][j + start_j];
            }
        }
        return ans;
    }
    //----------------------- Square matrix methods -----------------------
    T Trace() const requires(n == m) {
        T sum = T::ZERO;
        for (size_t i = 0; i < n && i < m; ++i) {
            sum += (*this)[i][i];
        }
    }

    Matrix Power(size_t indicator) const requires(n == m) {
        if (indicator == 0) return IdentityMatrix();
        if (indicator % 2 == 0) return (*this * *this).Power(indicator / 2);
        return *this * Power(indicator - 1);
    }

    Matrix Inverse() const requires(n == m) {
        auto ans = *this | IdentityMatrix();
        ans.Gauss();
        if (ans.template Slice<n, n>(0, 0) != IdentityMatrix()) {
            throw SingularMatrixException();
        }
        return ans.template Slice<n, n>(0, n);
    }

    Poly<T> CharPoly() const requires(n == m){
        Matrix<Fraction<Poly<T>>, n> ch;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                ch[i][j] -= Poly<T>((*this)[i][j]);
            }
        }
        for (size_t i = 0; i < n; ++i) {
            ch[i][i] += Poly<T>({T::ZERO, T::ONE});
        }
        return Det(ch).GetIntegerPart();
    }

    size_t rk() const {
        Matrix a = *this;
        a.Gauss();
        size_t j = 0;
        for (size_t i = 0; i < n; ++i) {
            while (j < m && a[i][j] == T::ZERO) {
                ++j;
            }
            if (j == m) {
               return i;
            }
        }
        return n;
    }

    static Matrix ScalarMatrix(const T& lambda) requires(n == m) {
        Matrix ans;
        for (size_t i = 0; i < n; ++i) {
            ans[i][i] = lambda;
        }
        return ans;
    }

    static Matrix IdentityMatrix() requires(n == m) {
        return ScalarMatrix(T::ONE);
    }
private:
    std::array<std::array<T, m>, n> data_;
};

template <RingWithOne T, size_t n, size_t m>
Matrix<T, n, m> operator*(const T& lambda, const Matrix<T, n, m>& matrix) {
    return matrix * lambda;
}

template <RingWithOne T, size_t n, size_t m>
std::ostream& operator<<(std::ostream& stream, const Matrix<T, n, m>& matrix) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            stream << matrix[i][j] << ' ';
        }
        stream << '\n';
    }
    return stream;
}



template <RingWithOne T, size_t n>
requires Field<T>
T Det(Matrix<T, n> matrix) {
    size_t start = 0;
    T ans = T::ONE;
    for (size_t j = 0; j < n; ++j) {
        size_t with_non_zero_coefficient = 0;
        bool found = false;
        for (size_t i = start; i < n; ++i) {
            if (matrix[i][j] != T::ZERO) {
                found = true;
                with_non_zero_coefficient = i;
                break;
            }
        }
        if (!found) {
            continue;
        }
        if (start != with_non_zero_coefficient) {
            ans = -ans;
        }
        swap(matrix[start], matrix[with_non_zero_coefficient]);
        for (size_t i = j + 1; i < n; ++i) {
            matrix[start][i] /= matrix[start][j];
        }
        ans *= matrix[start][j];
        matrix[start][j] = T::ONE;
        for (size_t i = 0; i < n; ++i) {
            if (i == start) {
                continue;
            }
            if (matrix[i][j] == T::ZERO) {
                continue;
            }
            auto lambda = -matrix[i][j];
            for (size_t k = j; k < n; ++k) {
                matrix[i][k] += lambda * matrix[start][k];
            }
        }
        ++start;
    }
    for (size_t i = 0; i < n; ++i) {
        ans *= matrix[i][i];
    }
    return ans;
}

template <RingWithOne T, size_t n>
requires (!Field<T>)
T Det(Matrix<T, n> matrix) {
    T ans = T::ZERO;
    for (const auto &indexes : AllPermutations(n)) {
        T cur = T::ONE;
        for (size_t i = 0; i < n; ++i) {
            cur *= matrix[i][indexes[i]];
        }
        if (indexes.sign() == -1) {
            ans -= cur;
        } else {
            ans += cur;
        }
    }
    return ans;
}

template <RingWithOne T, RingWithOne P, size_t n, size_t m>
Matrix<T, n, m> matrix_cast(const Matrix<P, n, m>& matrix) {
    Matrix<T, n, m> ans;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            ans[i][j] = static_cast<T>(matrix[i][j]);
        }
    }
    return ans;
}
