#pragma once

#include "exceptions.h"
#include "matrix.h"
#include "myconcepts.h"
#include "vector.h"

#include <algorithm>
#include <cstddef>
#include <vector>

template<typename T>
class VectorSpace {
public:
    VectorSpace(const std::vector<Vector<T>>& vectors) { // Linear span
        size_t m = vectors.size();
        size_t n = vectors[0].size();
        Matrix<T> a(n, m);
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < n; ++j) {
                a[j][i] = vectors[i][j];
            }
        }
        a.Gauss();
        size_t j = 0;
        for (size_t i = 0; i < n; ++i) {
            while (j < m && a[i][j] == T::ZERO()) {
                j++;
            }
            if (j == m) {
                break;
            }
            basis.push_back(vectors[j]);
        }
    }

    VectorSpace(Matrix<T> a) { // FSS
        size_t m = a.nsize();
        size_t n = a.msize();
        a.Gauss();
        std::vector<size_t> base_indexes;
        size_t j = 0;
        for (size_t i = 0; i < m; ++i) {
            while (j < n && a[i][j] == T::ZERO()) {
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
                Vector<T> new_vector(n);
                for (size_t j = 0; j < base_indexes.size(); ++j) {
                    new_vector[base_indexes[j]] = -a[j][i];
                }
                new_vector[i] = T::ONE();
                basis.push_back(new_vector);
            } else {
                ++current_j;
            }
        }
    }

    size_t dim() const {
        return basis.size();
    }

    size_t vector_size() const {
        return basis[0].size();
    }

    void MakeFullBasis() {
        size_t n = basis[0].size();
        Matrix<T> a = GetBasisMatrix();
        a.Gauss();
        size_t j = 0;
        size_t m = basis.size();
        for (size_t i = 0; i < n; ++i) {
            while (j < m && a[i][j] == T::ZERO()) {
                j++;
                if (j < m && a[i][j] == T::ZERO()) {
                    Vector<T> new_vector(n);
                    new_vector[j] = T::ONE();
                    basis.push_back(new_vector);
                }
            }
            if (j == m) {
                break;
            }
        }
        for (;j < n;++j) {
            Vector<T> new_vector(n);
            new_vector[j] = T::ONE();
            basis.push_back(new_vector);
        }
    }

    const std::vector<Vector<T>>& GetBasis() const {
        return basis;
    }

    std::vector<Vector<T>>& GetBasis() {
        return basis;
    }

    Matrix<T> GetBasisMatrix() const {
        size_t n = basis[0].size();
        size_t m = basis.size();
        Matrix<T> ans(n, m);
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < n; ++j) {
                ans[j][i] = basis[i][j];
            }
        }
        return ans;
    }
    Matrix<T> GetRowBasisMatrix() const {
        size_t n = basis[0].size();
        size_t m = basis.size();
        Matrix<T> ans(m, n);
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < n; ++j) {
                ans[i][j] = basis[i][j];
            }
        }
        return ans;
    }

    Matrix<T> GetEquasion() const {
        return VectorSpace(GetRowBasisMatrix()).GetRowBasisMatrix();
    }

    VectorSpace operator+(const VectorSpace& other) {
        if (vector_size() != other.vector_size()) {
            throw WrongSizeException();
        }
        size_t n1 = dim();
        size_t n2 = other.dim();
        std::vector<Vector<T>> sum_vectors(n1 + n2, Vector<T>(vector_size()));
        for (size_t i = 0; i < basis.size(); ++i) {
            sum_vectors[i] = basis[i];
        }
        for (size_t i = 0; i < other.basis.size(); ++i) {
            sum_vectors[i + n1] = other.basis[i];
        }
        return VectorSpace(sum_vectors);
    }

    VectorSpace inter(const VectorSpace& other) {
        if (vector_size() != other.vector_size()) {
            throw WrongSizeException();
        }
        Matrix<T> a = GetEquasion();
        Matrix<T> b = other.GetEquasion();
        Matrix<T> c(a.nsize() + b.nsize(), a.msize());
        for (size_t i = 0; i < a.nsize(); ++i) {
            for (size_t j = 0; j < a.msize(); ++j) {
                c[i][j] = a[i][j];
            }
        }
        for (size_t i = 0; i < b.nsize(); ++i) {
            for (size_t j = 0; j < b.msize(); ++j) {
                c[i + a.nsize()][j] = b[i][j];
            }
        }
        return VectorSpace(c);
    }

private:
    std::vector<Vector<T>> basis;
};
