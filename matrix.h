#pragma once

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
    Matrix() {}
    Matrix(const Matrix& other): data_(other.data_){}
    Matrix(size_t n, size_t m): data_(n, std::vector<T>(m, 0)) {}
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
    Matrix& operator=(const Matrix& other) {
        data_ = other.data_;
    }

    std::pair<size_t, size_t> size() const {
        return std::make_pair(data_.size(), data_[0].size());
    }

    Matrix operator+(const Matrix& other) const {
        if (other.size() != size()) throw WrongSizeException();
        auto [n, m] = size();
        Matrix ans(n, m);
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

    Matrix operator*(const Matrix& other) const {
        if (size().second != other.size().first) throw WrongSizeException();
        auto [n, m] = size();
        size_t k = other.size().second;
        Matrix a(n, k);
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
                if (data_[i][j] != 0) {
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
            data_[start][j] = 1;
            for (size_t i = 0; i < n; ++i) {
                if (i == start) {
                    continue;
                }
                if (data_[i][j] == 0) {
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
    }

    T Trace() const {
        auto [n, m] = size();
        T sum = 0;
        for (int i = 0; i < n && i < m; ++i) {
            sum += data_[i][i];
        }
    }

protected:
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
    SquareMatrix(int n): Matrix<T>(n, n) {}
    int size() const {
        return size().first;
    }

    SquareMatrix Power(size_t indicator) const {
        if (indicator == 0) return IdentityMatrix();
        if (indicator % 2 == 0) return (*this * *this).Power(indicator / 2);
        return Power(indicator - 1) * *this;
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

    static SquareMatrix ScalarMatrix(int n, const T& lambda) {
        SquareMatrix ans(n);
        for (int i = 0; i < n; ++i) {
            ans[i][i] = lambda;
        }
        return ans;
    }

    static SquareMatrix IdentityMatrix(int n) {
        return ScalarMatrix(n, 1);
    }
};
