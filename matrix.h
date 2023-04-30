#pragma once

#include "fraction.h"
#include "exceptions.h"
#include "permutation.h"
#include "poly.h"
#include "myconcepts.h"

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>

template <RingWithOne T>
class Matrix {
public:
    Matrix() = delete;
    Matrix(size_t n, size_t m) : data_(n, std::vector<T>(m, T::ZERO())) {
    }
    Matrix(std::pair<size_t, size_t> size) : Matrix(size.first, size.second) {
    }
    Matrix(const Matrix& other) = default;
    Matrix(Matrix&& other) = default;
    Matrix& operator=(const Matrix& other) = default;
    Matrix& operator=(Matrix&& other) = default;

    Matrix(std::initializer_list<std::initializer_list<T>> data) {
        data_.assign(data.size(), std::vector<T>(data.begin() -> size()));
        size_t i = 0;
        for (const auto& row : data) {
            if (row.size() != data.begin() -> size()) {
                throw WrongSizeException();
            }
            size_t j = 0;
            for (const auto& elem : row) {
                data_[i][j] = elem;
                ++j;
            }
            ++i;
        }
    }

    bool operator==(const Matrix& other) const = default;
    bool operator!=(const Matrix& other) const = default;

    size_t nsize() const {
        return data_.size();
    }
    size_t msize() const {
        return data_[0].size();
    }
    std::pair<size_t, size_t> size() const {
        return std::make_pair(nsize(), msize());
    }

    Matrix operator+(const Matrix& other) const {
        if (other.size() != size()) {
            throw WrongSizeException();
        }
        size_t n = nsize();
        size_t m = msize();
        Matrix ans(size());
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                ans[i][j] = data_[i][j] + other[i][j];
            }
        }
        return ans;
    }

    Matrix& operator+=(const Matrix& other) {
        return *this = *this + other;
    }

    Matrix operator-() const {
        Matrix ans = *this;
        size_t n = nsize();
        size_t m = msize();
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

    Matrix<T> operator*(const Matrix<T>& other) const {
        if (msize() != other.nsize()) {
            throw WrongSizeException();
        }
        size_t n = nsize();
        size_t m = msize();
        size_t k = other.msize();
        Matrix<T> a(n, k);
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
        Matrix ans(size());
        size_t n = nsize();
        size_t m = msize();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                ans[i][j] = lambda * data_[i][j];
            }
        }
        return ans;
    }

    Matrix& operator*=(const T& lambda) {
        size_t n = nsize();
        size_t m = msize();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                data_[i][j] *= lambda;
            }
        }
        return *this;
    }

    // writes other to the right from *this
    Matrix<T> operator|(const Matrix<T>& other) const {
        if (nsize() != other.nsize()) {
            throw WrongSizeException();
        }
        size_t n = nsize();
        size_t m = msize();
        size_t k = other.msize();
        Matrix<T> a(n, m + k);
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

    std::vector<T>& operator[](size_t pos) {
        return data_[pos];
    }
    const std::vector<T>& operator[](size_t pos) const {
        return data_[pos];
    }

    void Gauss() {
        size_t n = nsize();
        size_t m = msize();
        size_t start = 0;
        for (size_t j = 0; j < m; ++j) {
            size_t with_non_zero_coefficient = 0;
            bool found = false;
            for (size_t i = start; i < n; ++i) {
                if (data_[i][j] != T::ZERO()) {
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
            data_[start][j] = T::ONE();
            for (size_t i = 0; i < n; ++i) {
                if (i == start) {
                    continue;
                }
                if (data_[i][j] == T::ZERO()) {
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
        size_t n = nsize();
        size_t m = msize();
        Matrix<T> ans(m, n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                ans[j][i] = data_[i][j];
            }
        }
        return ans;
    }

    Matrix<T> Slice(size_t len_i, size_t len_j, size_t start_i = 0, size_t start_j = 0) const {
        size_t n = nsize();
        size_t m = msize();
        if (len_i + start_i > n || len_j + start_j > m)
            throw WrongSizeException();
        Matrix<T> ans(len_i, len_j);
        for (size_t i = 0; i < len_i; ++i) {
            for (size_t j = 0; j < len_j; ++j) {
                ans[i][j] = data_[i + start_i][j + start_j];
            }
        }
        return ans;
    }
    //----------------------- Square matrix methods -----------------------
    T Trace() const {
        size_t n = nsize();
        size_t m = msize();
        if (n != m) {
            throw WrongSizeException();
        }
        T sum = T::ZERO();
        for (size_t i = 0; i < n && i < m; ++i) {
            sum += (*this)[i][i];
        }
    }

    Matrix Power(size_t indicator) const {
        if (nsize() != msize()) {
            throw WrongSizeException();
        }
        if (indicator == 0)
            return IdentityMatrix(nsize());
        if (indicator % 2 == 0)
            return (*this * *this).Power(indicator / 2);
        return *this * Power(indicator - 1);
    }

    Matrix Inverse() const {
        size_t n = nsize();
        size_t m = msize();
        if (n != m) {
            throw WrongSizeException();
        }
        auto ans = *this | IdentityMatrix(n);
        ans.Gauss();
        if (ans.Slice(n, n, 0, 0) != IdentityMatrix(n)) {
            throw SingularMatrixException();
        }
        return ans.Slice(n, n, 0, n);
    }

    Poly<T> CharPoly() const requires(Field<T>) {
        size_t n = nsize();
        size_t m = msize();
        if (n != m) {
            throw WrongSizeException();
        }
        Matrix<Fraction<Poly<T>>> ch(n, n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                ch[i][j] -= Poly<T>((*this)[i][j]);
            }
        }
        for (size_t i = 0; i < n; ++i) {
            ch[i][i] += Poly<T>({T::ZERO(), T::ONE()});
        }
        return Det(ch).GetIntegerPart();
    }

    Poly<T> CharPoly() const requires(!Field<T>) {
        size_t n = nsize();
        size_t m = msize();
        if (n != m) {
            throw WrongSizeException();
        }
        Matrix<Poly<T>> ch(n, n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                ch[i][j] -= Poly<T>((*this)[i][j]);
            }
        }
        for (size_t i = 0; i < n; ++i) {
            ch[i][i] += Poly<T>({T::ZERO(), T::ONE()});
        }
        return Det(ch);
    }

    size_t rk() const {
        size_t n = nsize();
        size_t m = msize();
        if (n != m) {
            throw WrongSizeException();
        }
        Matrix a = *this;
        a.Gauss();
        size_t j = 0;
        for (size_t i = 0; i < n; ++i) {
            while (j < m && a[i][j] == T::ZERO()) {
                ++j;
            }
            if (j == m) {
                return i;
            }
        }
        return n;
    }

    static Matrix ScalarMatrix(size_t n, const T& lambda) {
        Matrix ans(n, n);
        for (size_t i = 0; i < n; ++i) {
            ans[i][i] = lambda;
        }
        return ans;
    }

    static Matrix IdentityMatrix(size_t n) {
        return ScalarMatrix(n, T::ONE());
    }

private:
    std::vector<std::vector<T>> data_;
};

template <RingWithOne T>
Matrix<T> operator*(const T& lambda, const Matrix<T>& matrix) {
    return matrix * lambda;
}

template <RingWithOne T>
std::ostream& operator<<(std::ostream& stream, const Matrix<T>& matrix) {
    std::vector<size_t> maxsizes(matrix.msize());
    for (size_t i = 0; i < matrix.nsize(); ++i) {
        for (size_t j = 0; j < matrix.msize(); ++j) {
            std::stringstream tmp;
            tmp << matrix[i][j];
            maxsizes[j] = std::max(maxsizes[j], tmp.str().size());
            // stream << streams[i][j].str() << ' ';
            // stream << matrix[i][j] << ' ';
        }
        // stream << '\n';
    }
    for (size_t i = 0; i < matrix.nsize(); ++i) {
        if (i == 0) {
            stream << "⎛ ";
            // stream << "⎡ ";
        } else if (i == matrix.nsize() - 1) {
            stream << "⎝ ";
            // stream << "⎣ ";
        } else {
            stream << "⎜ ";
            // stream << "⎢ ";
        }
        for (size_t j = 0; j < matrix.msize(); ++j) {
            std::stringstream tmp;
            tmp << matrix[i][j];
            stream << (j ? " " : "") << std::setw(static_cast<int>(maxsizes[j])) << tmp.str();
        }
        if (i == 0) {
            stream << " ⎞\n";
            // stream << " ⎤\n";
        } else if (i == matrix.nsize() - 1) {
            stream << " ⎠\n";
            // stream << " ⎦\n";
        } else {
            stream << " ⎟\n";
            // stream << " ⎥\n";
        }
    }
    return stream;
}

template <RingWithOne T>
    requires Field<T>
T Det(Matrix<T> matrix) {
    size_t n = matrix.nsize();
    size_t m = matrix.msize();
    if (n != m) {
        throw WrongSizeException();
    }
    size_t start = 0;
    T ans = T::ONE();
    for (size_t j = 0; j < n; ++j) {
        size_t with_non_zero_coefficient = 0;
        bool found = false;
        for (size_t i = start; i < n; ++i) {
            if (matrix[i][j] != T::ZERO()) {
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
        matrix[start][j] = T::ONE();
        for (size_t i = 0; i < n; ++i) {
            if (i == start) {
                continue;
            }
            if (matrix[i][j] == T::ZERO()) {
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

template <RingWithOne T>
    requires(!Field<T>)
T Det(Matrix<T> matrix) {
    size_t n = matrix.nsize();
    size_t m = matrix.msize();
    if (n != m) {
        throw WrongSizeException();
    }
    T ans = T::ZERO();
    for (const auto& indexes : AllPermutations(n)) {
        T cur = T::ONE();
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

template <RingWithOne T, RingWithOne P>
Matrix<T> matrix_cast(const Matrix<P>& matrix) {
    size_t n = matrix.nsize();
    size_t m = matrix.msize();
    Matrix<T> ans(matrix.size());
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            ans[i][j] = static_cast<T>(matrix[i][j]);
        }
    }
    return ans;
}
