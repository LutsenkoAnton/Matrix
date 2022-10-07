#include "matrix.h"

template<size_t n, size_t m>
Matrix<n, m>::Matrix() {
}

template<size_t n, size_t m>
Matrix<n, m>::Matrix(const Matrix& other) : data_(other.data_) {
}

template<size_t n, size_t m>
Matrix<n, m>::Matrix(std::array<std::array<Fraction, m>, n> data) : data_(data) {
}

template<size_t n, size_t m>
Matrix<n, m>::Matrix(std::initializer_list<std::initializer_list<Fraction>> data) {
    size_t i = 0;
    for (const auto& row : data) {
        size_t j = 0;
        for (const auto& elem : row) {
            data_[i][j] = elem;
            ++j;
        }
        ++i;
    }
}

template<size_t n, size_t m>
Matrix<n, m>& Matrix<n, m>::operator=(const Matrix& other) {
    data_ = other.data_;
    return *this;
}

template<size_t n, size_t m>
Matrix<n, m> Matrix<n, m>::operator+(const Matrix& other) const {
    Matrix a;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            a[i][j] = data_[i][j] + other[i][j];
        }
    }
}

template<size_t n, size_t m>
Matrix<n, m>& Matrix<n, m>::operator+=(const Matrix& other) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            data_[i][j] += other[i][j];
        }
    }
    return *this;
}

template<size_t n, size_t m>
template<size_t k>
Matrix<n, k> Matrix<n, m>::operator*(const Matrix<m, k>& other) const {
    Matrix<n, k> a;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            for (size_t l = 0; l < k; ++l) {
                a[i][l] += data_[i][j] * other[j][l];
            }
        }
    }
    return a;
}

template<size_t n, size_t m>
Matrix<n, m> Matrix<n, m>::operator*(const Fraction lambda) const {
    Matrix a(*this);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            a[i][j] *= lambda;
        }
    }
    return a;
}

template<size_t n, size_t m>
Matrix<n, m>& Matrix<n, m>::operator*=(const Fraction lambda) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            data_[i][j] *= lambda;
        }
    }
    return *this;
}

template<size_t n, size_t m>
template<size_t k>
Matrix<n, m + k> Matrix<n, m>::operator|(const Matrix<n, k>& other) const {
    Matrix<n, m + k> a;
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

template<size_t n, size_t m>
std::array<Fraction, m>& Matrix<n, m>::operator[](size_t pos) {
    return data_[pos];
}

template<size_t n, size_t m>
std::array<Fraction, m> Matrix<n, m>::operator[](size_t pos) const {
    return data_[pos];
}

template<size_t n, size_t m>
void Matrix<n, m>::Gauss() {
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

template<size_t n, size_t m>
Matrix<m, n> Matrix<n, m>::Transpose() const {
    Matrix<m, n> a;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            a[j][i] = data_[i][j];
        }
    }
    return a;
}

template<size_t n, size_t m>
template<size_t x, size_t y>
void Matrix<n, m>::Slice(Matrix<x, y>& to, size_t start_i, size_t start_j) const {
    for (size_t i = 0; i < x; ++i) {
        for (size_t j = 0; j < y; ++j) {
            to[i][j] = data_[i + start_i][j + start_j];
        }
    }
}

template<size_t n, size_t m>
Fraction Matrix<n, m>::Trace() const {
    Fraction ans;
    for (size_t i = 0; i < std::min(n, m); ++i) {
        ans += data_[i][i];
    }
    return ans;
}

template <size_t n, size_t m>
Matrix<n, m> operator*(const Fraction lambda, const Matrix<n, m>& matrix) {
    return matrix * lambda;
}

template <size_t n, size_t m>
std::ostream& operator<<(std::ostream& stream, const Matrix<n, m>& matrix) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            stream << matrix[i][j] << ' ';
        }
        stream << '\n';
    }
    return stream;
}

template <size_t n>
SquareMatrix<n>::SquareMatrix() {
}
template <size_t n>
SquareMatrix<n>::SquareMatrix(Fraction lambda) : Matrix<n, n>(lambda * IdentityMatrix()) {
}

template <size_t n>
SquareMatrix<n> SquareMatrix<n>::Power(size_t indicator) const {
    SquareMatrix a = IdentityMatrix();
    for (size_t i = 0; i < indicator; ++i) {
        a *= (*this);
    }
    return a;
}

template <size_t n>
SquareMatrix<n> SquareMatrix<n>::Inverse() const {
    auto a = *this | IdentityMatrix();
    a.Gauss();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i != j && a[i][j] != 0) {
                return Fraction(0);
            }
            if (i == j && a[i][i] != 1) {
                return Fraction(0);
            }
        }
    }
    SquareMatrix ans;
    this->Slice(a, ans, 0, n);
    return ans;
}

template <size_t n>
SquareMatrix<n> SquareMatrix<n>::IdentityMatrix() {
    SquareMatrix a;
    for (size_t i = 0; i < n; ++i) {
        a[i][i] = 1;
    }
    return a;
}
