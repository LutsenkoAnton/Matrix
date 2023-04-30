#pragma once

#include "matrix.h"
#include "myconcepts.h"
#include "poly.h"
#include "vector.h"
#include "vector_space.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <exception>
#include <iterator>
#include <utility>
#include <vector>

template <Field T>
class LinearOperator {
public:
    LinearOperator(size_t n): data_(n, n) {}
    LinearOperator(const Matrix<T>& data) : data_(data) {
    }
    LinearOperator(const VectorSpace<T>& im, const VectorSpace<T>& ker)
        : LinearOperator(im.GetBasisMatrix() * ker.GetEquasion()) {
    }
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

    VectorSpace<T> Im() const {
        size_t n = data_.nsize();
        std::vector<Vector<T>> vectors(n, Vector<T>(n));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                vectors[i][j] = data_[j][i];
            }
        }
        return VectorSpace<T>(vectors);
    }

    VectorSpace<T> ker() const {
        return VectorSpace<T>(data_);
    }

    Matrix<T> GetData() const {
        return data_;
    }

    std::vector<std::pair<T, size_t>> GetJNFBlocks() const {
        std::vector<std::pair<T, size_t>> jordan_blocks;
        auto solutions = data_.CharPoly().try_solve();
        auto e = Matrix<T>::IdentityMatrix(data_.nsize());
        for (auto [lambda, n_i] : solutions) {
            size_t k = 1;
            while (n_i > 0) {
                size_t cnt_k = (data_ - lambda * e).Power(k + 1).rk() + (data_ - lambda * e).Power(k - 1).rk() - 2 * (data_ - lambda * e).Power(k).rk();
                for (size_t i = 0; i < cnt_k; ++i) {
                    jordan_blocks.emplace_back(lambda, k);
                }
                n_i -= cnt_k * k;
                ++k;
            }
            assert(n_i == 0);
        }
        return jordan_blocks;
    }

    Matrix<T> GetJNF() const {
        Matrix<T> ans(data_.size());
        size_t offset = 0;
        auto jordan_blocks = GetJNFBlocks();
        for (auto [lambda, k] : jordan_blocks) {
            for (size_t i = 0; i < k; ++i) {
                ans[i + offset][i + offset] = lambda;
            }
            for (size_t i = 1; i < k; ++i) {
                ans[i + offset - 1][i + offset] = T::ONE();
            }
            offset += k;
        }
        return ans;
    }

    Matrix<T> GetJordanBasis() const {
        auto jordan_blocks = GetJNFBlocks();
        std::reverse(jordan_blocks.begin(), jordan_blocks.end());
        size_t i = 0;
        auto e = Matrix<T>::IdentityMatrix(data_.nsize());
        std::vector<Vector<T>> basis;
        while (i < jordan_blocks.size()) {
            auto [lambda, n] = jordan_blocks[i];
            std::vector<std::vector<Vector<T>>> chains;
            std::vector<Vector<T>> current_vectors;
            for (size_t j = n; j > 0; --j) {
                for (auto& elem : current_vectors) {
                    elem = (data_ - lambda * e) * elem;
                }
                size_t current_size = current_vectors.size();
                auto kerj = LinearOperator<T>((data_ - lambda * e).Power(j)).ker();
                auto kerjm1 = LinearOperator<T>((data_ - lambda * e).Power(j - 1)).ker();
                std::vector<Vector<T>> to_filter;
                for (const auto &elem : kerjm1.GetBasis()) {
                    to_filter.push_back(elem);
                }
                for (const auto &elem : kerj.GetBasis()) {
                    to_filter.push_back(elem);
                }
                auto filtered = VectorSpace<T>(to_filter).GetBasis();
                for (size_t j2 = kerjm1.dim(); j2 < filtered.size(); ++j2) {
                    current_vectors.push_back(filtered[j2]);
                }
                current_vectors = VectorSpace<T>(current_vectors).GetBasis();
                for (ssize_t j2 = current_vectors.size() - 1; j2 >= static_cast<ssize_t>(current_size); --j2) {
                    std::vector<Vector<T>> new_chain(j, Vector<T>(data_.nsize()));
                    new_chain[0] = current_vectors[j2];
                    for (size_t k = 1; k < j; ++k) {
                        new_chain[k] = (data_ - lambda * e) * new_chain[k - 1];
                    }
                    // reverse(new_chain.begin(), new_chain.end());
                    chains.push_back(new_chain);
                }
            }
            // std::reverse(chains.begin(), chains.end());
            for (const auto& chain : chains) {
                for (const auto& vec : chain) {
                    basis.push_back(vec);
                }
            }
            while (jordan_blocks.size() > i && lambda == jordan_blocks[i].first) {++i;}
        }
        std::reverse(basis.begin(), basis.end());
        return VectorSpace(basis).GetBasisMatrix();
    }

    static LinearOperator ZERO(size_t n) {
        return LinearOperator(n);
    }
    static LinearOperator ONE(size_t n) {
        return LinearOperator(Matrix<T>::IdentityMatrix(n));
    }

private:
    Matrix<T> data_;
};
