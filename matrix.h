#pragma once

#include "linear_combination.h"
#include "permutation.h"
#include "poly.h"
#include "rational_function.h"
#include "variable.h"

#include <cassert>
#include <initializer_list>
#include <iostream>
#include <exception>
#include <numeric>
#include <vector>

class WrongSizeException: public std::exception {
public:
    const char* what() const noexcept override {
        return "Wrong sizes of matrixes";
    }
};

template <typename T>
class Matrix {
public:
    Matrix() = delete;
    Matrix(const Matrix& other) = default;
    Matrix(Matrix&& other) = default;
    Matrix& operator=(const Matrix& other) = default;
    Matrix& operator=(Matrix&& other) = default;

    Matrix(size_t n, size_t m): data_(n, std::vector<T>(m, static_cast<T>(0))) {}
    explicit Matrix(std::vector<std::vector<T>> data): data_(data) {
        for (int i = 0; i < data_.size(); ++i) {
            if (data_[i].size() != data_[0].size()) throw WrongSizeException();
        }
    }
    Matrix(std::initializer_list<std::initializer_list<T>> data) {
        for (const auto& row : data) {
            data_.emplace_back();
            for (const auto& elem : row) {
                data_.back().push_back(elem);
            }
        }
        for (int i = 0; i < data_.size(); ++i) {
            if (data_[i].size() != data_[0].size()) throw WrongSizeException();
        }
    }

    bool operator==(const Matrix& other) const = default;
    bool operator!=(const Matrix& other) const = default;

    std::pair<size_t, size_t> size() const {
        return std::make_pair(data_.size(), data_[0].size());
    }

    //decltype for adding matrix of LinearCombinations and matrix of coefficients
    template<typename P>
    Matrix<decltype(T() + P())>  operator+(const Matrix<P>& other) const {
        if (other.size() != size()) throw WrongSizeException();
        auto [n, m] = size();
        Matrix<decltype(T() + P())>  ans(n, m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                ans[i][j] = data_[i][j] + other[i][j];
            }
        }
        return ans;
    }

    Matrix& operator+=(const Matrix& other) {
        if (other.size() != size()) throw WrongSizeException();
        auto [n, m] = size();
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                data_[i][j] += other[i][j];
            }
        }
        return *this;
    }

    Matrix operator-() const {
        auto [n, m] = size();
        Matrix ans(data_);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                ans[i][j] = -ans[i][j];
            }
        }
        return ans;
    }

    //decltype for substracting matrix of LinearCombinations from matrix of coefficients and vice versa
    template<typename P>
    Matrix<decltype(T() - P())> operator-(const Matrix<P>& other) const {
        return *this + -other;
    }

    Matrix& operator-=(const Matrix& other) const {
        return *this = *this - other;
    }

    //decltype for multiplying matrix of LinearCombinations by matrix of coefficients
    template<typename P>
    Matrix<decltype(T() * P())> operator*(const Matrix<P>& other) const {
        if (size().second != other.size().first) throw WrongSizeException();
        auto [n, m] = size();
        size_t k = other.size().second;
        Matrix<decltype(T() * P())> a(n, k);
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
        auto [n, m] = size();
        Matrix ans(n, m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                ans[i][j] = lambda * data_[i][j];
            }
        }
        return ans;
    }

    Matrix& operator*=(const T& lambda) {
        auto [n, m] = size();
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                data_[i][j] *= lambda;
            }
        }
    }

    // writes other to the right from *this
    Matrix operator|(const Matrix& other) const {
        if (size().first != other.size().first) throw WrongSizeException();
        auto [n, m] = size();
        size_t k = other.size().second;
        Matrix a(n, m + k);
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
    std::vector<T> operator[](size_t pos) const {
        return data_[pos];
    }

    void Gauss() {
        auto [n, m] = size();
        size_t start = 0;
        for (size_t j = 0; j < m; ++j) {
            size_t with_non_zero_coefficient = -1;
            bool found = false;
            for (size_t i = start; i < n; ++i) {
                if (data_[i][j] != T(0)) {
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
            data_[start][j] = T(1);
            for (size_t i = 0; i < n; ++i) {
                if (i == start) {
                    continue;
                }
                if (data_[i][j] == T(0)) {
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
        auto [n, m] = size();
        Matrix ans(m, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                ans[j][i] = data_[i][j];
            }
        }
        return ans;
    }

    Matrix Slice(size_t len_i, size_t len_j, size_t start_i = 0, size_t start_j = 0) const {
        auto [n, m] = size();
        if (len_i + start_i > n || len_j + start_j > m) throw WrongSizeException();
        Matrix ans(len_i, len_j);
        for (int i = 0; i < len_i; ++i) {
            for (int j = 0; j < len_j; ++j) {
                ans[i][j] = data_[i + start_i][j + start_j];
            }
        }
        return ans;
    }

    T Trace() const {
        auto [n, m] = size();
        T sum = 0;
        for (int i = 0; i < n && i < m; ++i) {
            sum += data_[i][i];
        }
    }

    static Matrix<LinearCombination<T>> MakeVariable(size_t n, size_t m) {
        Matrix<LinearCombination<T>> res(n, m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                res[i][j] = Variable();
            }
        }
        return res;
    }

private:
    std::vector<std::vector<T>> data_;
};

template <typename T>
Matrix<T> operator*(const T& lambda, const Matrix<T>& matrix) {
    return matrix * lambda;
}

template <typename T>
std::ostream& operator<<(std::ostream& stream, const Matrix<T>& matrix) {
    auto [n, m] = matrix.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            stream << matrix[i][j] << ' ';
        }
        stream << '\n';
    }
    return stream;
}

class SingularMatrixException: public std::exception {
public:
    const char* what() const noexcept override {
        return "Matrix is singular";
    }
};

template <typename T>
class SquareMatrix : public Matrix<T> {
public:
    SquareMatrix() = delete;
    SquareMatrix(const SquareMatrix& other) = default;
    SquareMatrix(SquareMatrix&& other) = default;
    SquareMatrix& operator=(const SquareMatrix& other) = default;
    SquareMatrix& operator=(SquareMatrix&& other) = default;
    SquareMatrix(int n): Matrix<T>(n, n) {}
    explicit SquareMatrix(std::vector<std::vector<T>> data): Matrix<T>(data) {}
    SquareMatrix(std::initializer_list<std::initializer_list<T>> data): Matrix<T>(data) {}
    SquareMatrix(const Matrix<T>& m): Matrix<T>(m) {
        if (m.size().first != m.size().second) throw WrongSizeException();
    }
    
    int size() const {
        return Matrix<T>::size().first;
    }
    template<typename P>
    SquareMatrix<decltype(T() + P())> operator+(const SquareMatrix<P>& other) const {
        return Matrix<T>::operator+(other);
    }
    template<typename P>
    SquareMatrix<decltype(T() - P())> operator-(const SquareMatrix<P>& other) const {
        return Matrix<T>::operator-(other);
    }
    template<typename P>
    SquareMatrix<decltype(T() * P())> operator*(const SquareMatrix<P>& other) const {
        return Matrix<T>::operator*(other);
    }
    SquareMatrix operator-() const {
        return Matrix<T>::operator-();
    }
    template<typename P>
    Matrix<decltype(T() + P())> operator+(const Matrix<P>& other) const {
        return Matrix<T>::operator+(other);
    }
    template<typename P>
    Matrix<decltype(T() - P())> operator-(const Matrix<P>& other) const {
        return Matrix<T>::operator-(other);
    }
    template<typename P>
    Matrix<decltype(T() * P())> operator*(const Matrix<P>& other) const {
        return Matrix<T>::operator*(other);
    }

    SquareMatrix Power(size_t indicator) const {
        if (indicator == 0) return IdentityMatrix(size());
        if (indicator % 2 == 0) return (*this * *this).Power(indicator / 2);
        return *this * Power(indicator - 1);
    }

    SquareMatrix Inverse() const {
        size_t n = size();
        auto ans = *this | IdentityMatrix(n);
        ans.Gauss();
        if (ans.Slice(n, n) != IdentityMatrix(n)) {
            throw SingularMatrixException();
        }
        return ans.Slice(n, n, 0, n);
    }

    Poly<T> CharPoly() const {
        int n = size();
        SquareMatrix<RationalFunction<T>> ch(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                ch[i][j] -= Poly<T>((*this)[i][j]);
            }
        }
        for (size_t i = 0; i < n; ++i) {
            ch[i][i] += Poly<T>({T(0), T(1)});
        }
        return Det(ch).force_devide();
    }

    Poly<T> PermCharPoly() const {
        int n = size();
        SquareMatrix<Poly<T>> ch(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                ch[i][j] -= (*this)[i][j];
            }
        }
        for (size_t i = 0; i < n; ++i) {
            ch[i][i] += Poly<T>({T(0), T(1)});
        }
        return PermDet(ch);
    }

    static SquareMatrix ScalarMatrix(int n, const T& lambda) {
        SquareMatrix ans(n);
        for (int i = 0; i < n; ++i) {
            ans[i][i] = lambda;
        }
        return ans;
    }

    static SquareMatrix IdentityMatrix(int n) {
        return ScalarMatrix(n, T(1));
    }
};

template <typename T>
T Det(SquareMatrix<T> matrix) {
    size_t n = matrix.size();
    size_t start = 0;
    T ans = T(1);
    for (size_t j = 0; j < n; ++j) {
        size_t with_non_zero_coefficient = 0;
        bool found = false;
        for (size_t i = start; i < n; ++i) {
            if (matrix[i][j] != T(0)) {
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
        matrix[start][j] = T(1);
        for (size_t i = 0; i < n; ++i) {
            if (i == start) {
                continue;
            }
            if (matrix[i][j] == T(0)) {
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

template <typename T>
T PermDet(SquareMatrix<T> matrix) {
    size_t n = matrix.size();
    T ans;
    for (const auto &indexes : AllPermutations(n)) {
        size_t cnt = 0;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < i; ++j) {
                cnt += (indexes[j] > indexes[i]);
            }
        }
        T cur = T(1);
        for (size_t i = 0; i < n; ++i) {
            cur *= matrix[i][indexes[i]];
        }
        if (cnt % 2 == 0) {
            ans += cur;
        } else {
            ans -= cur;
        }
    }
    return ans;
}

class NotSolvableException: std::exception {
public:
    const char* what() const noexcept override {
        return "The system cannot be solved";
    }
};

template <typename T>
std::vector<LinearCombination<T>> Solve(const std::vector<LinearCombination<T>>& combinations) {
    size_t n = Variable::count;
    Matrix<T> slv(combinations.size(), n + 1);
    for (size_t i = 0; i < combinations.size(); ++i) {
        for (const auto& [var, coeff] : combinations[i]) {
            if (var.id == 0) {
                slv[i][n] = coeff;
            } else {
                slv[i][var.id - 1] = coeff;
            }
        }
    }
    slv.Gauss();
    std::vector<LinearCombination<T>> ans(n);
    for (size_t i = 0; i < n; ++i) {
        ans[i] = Variable(i + 1);
    }
    for (size_t i = 0; i < combinations.size(); ++i) {
        size_t j = 0;
        while (j <= n && slv[i][j] == 0) {
            j++;
        }
        if (j == n + 1) break;
        if (j == n) throw NotSolvableException();
        ans[j] -= Variable(j + 1);
        for (size_t k = j + 1; k < n; ++k) {
            ans[j] -= LinearCombination<T>(Variable(k + 1)) * slv[i][k];
        }
        ans[j] -= slv[i][n];
    }
    return ans;
}

template <typename T>
std::vector<LinearCombination<T>> Solve(const Matrix<LinearCombination<T>>& A) {
    std::vector<LinearCombination<T>> elements;
    auto [n, m] = A.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            elements.push_back(A[i][j]);
        }
    }
    return Solve(elements);
}

template <typename T, typename P>
Matrix<T> matrix_cast(const Matrix<P>& matrix) {
    auto [n, m] = matrix.size();
    Matrix<T> ans(matrix.size());
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            ans[i][j] = static_cast<T>(matrix[i][j]);
        }
    }
    return ans;
}

template <typename T, typename P>
SquareMatrix<T> matrix_cast(const SquareMatrix<P>& matrix) {
    auto n = matrix.size();
    SquareMatrix<T> ans(matrix.size());
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            ans[i][j] = static_cast<T>(matrix[i][j]);
        }
    }
    return ans;
}
