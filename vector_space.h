#pragma once

#include "matrix.h"
#include "myconcepts.h"
#include "vector.h"

#include <algorithm>
#include <cstddef>
#include <vector>

template<Field T, size_t n>
class VectorSpace {
public:
    template<size_t m>
    VectorSpace(const std::array<Vector<T, n>, m>& vectors) { // Linear span
        Matrix<T, n, m> a;
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < n; ++j) {
                a[j][i] = vectors[i][j];
            }
        }
        a.Gauss();
        size_t j = 0;
        for (size_t i = 0; i < n; ++i) {
            while (j < m && a[i][j] == T::ZERO) {
                j++;
            }
            if (j == m) {
                break;
            }
            basis.push_back(vectors[j]);
        }
    }

    template<size_t m>
    VectorSpace(Matrix<T, m, n> a) { // FSS
        a.Gauss();
        std::vector<size_t> base_indexes;
        size_t j = 0;
        for (size_t i = 0; i < m; ++i) {
            while (j < n && a[i][j] == T::ZERO) {
                j++;
            }
            if (j == n) {
                break;
            }
            base_indexes.push_back(j);
        }
        size_t current_j = 0;
        for (size_t i = 0; i < n; ++i) {
            if (current_j == base_indexes.size() || i != base_indexes[current_j]) {
                Vector<T, n> new_vector;
                for (size_t j = 0; j < base_indexes.size(); ++j) {
                    new_vector[base_indexes[j]] = -a[j][i];
                }
                new_vector[i] = T::ONE;
                basis.push_back(new_vector);
            } else {
                ++current_j;
            }
        }
    }

    size_t dim() const {
        return basis.size();
    }

    void MakeFullBasis() {
        Matrix<T, n> a = GetBasisMatrix();
        a.Gauss();
        size_t j = 0;
        for (size_t i = 0; i < n; ++i) {
            while (j < n && a[i][j] == 0) {
                j++;
                if (a[i][j] == 0) {
                    Vector<T, n> new_vector;
                    new_vector[j] = T::ONE;
                    basis.push_back(new_vector);
                }
            }
            if (j == n) {
                break;
            }
        }
    }

    const std::vector<Vector<T, n>>& GetBasis() const {
        return basis;
    }

    std::vector<Vector<T, n>>& GetBasis() {
        return basis;
    }

    Matrix<T, n> GetBasisMatrix() const {
        Matrix<T, n> ans;
        for (size_t i = 0; i < basis.size(); ++i) {
            for (size_t j = 0; j < n; ++j) {
                ans[j][i] = basis[i][j];
            }
        }
        return ans;
    }
    Matrix<T, n> GetRowBasisMatrix() const {
        Matrix<T, n> ans;
        for (size_t i = 0; i < basis.size(); ++i) {
            for (size_t j = 0; j < n; ++j) {
                ans[i][j] = basis[i][j];
            }
        }
        return ans;
    }

    Matrix<T, n> GetEquasion() const {
        return VectorSpace(GetRowBasisMatrix()).GetRowBasisMatrix();
    }

    VectorSpace operator+(const VectorSpace& other) {
        std::array<Vector<T, n>, 2 * n> sum_vectors;
        for (size_t i = 0; i < basis.size(); ++i) {
            sum_vectors[i] = basis[i];
        }
        for (size_t i = 0; i < other.basis.size(); ++i) {
            sum_vectors[i + n] = other.basis[i];
        }
        return VectorSpace(sum_vectors);
    }

    VectorSpace inter(const VectorSpace& other) {
        Matrix<T, n> a = GetEquasion();
        Matrix<T, n> b = other.GetEquasion();
        Matrix<T, 2 * n, n> c;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                c[i][j] = a[i][j];
            }
        }
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                c[i + n][j] = b[i][j];
            }
        }
        return VectorSpace(c);
    }

private:
    std::vector<Vector<T, n>> basis;
};
