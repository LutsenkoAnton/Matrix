#pragma once

#include <algorithm>
#include <array>
#include <initializer_list>
#include <iostream>
#include <numeric>
#include <vector>

class Permutation {
public:
    Permutation(size_t n) {
        data_.resize(n);
        std::iota(data_.begin(), data_.end(), 0);
    }
    Permutation(const std::vector<size_t>& data) : data_(data) {}
    Permutation(const std::initializer_list<size_t>& data) {
        *this = data;
    }
    Permutation(const std::initializer_list<std::initializer_list<size_t>>& cycles) {
        *this = cycles;
    }
    Permutation(const Permutation& other) = default;
    Permutation(Permutation&& other) = default;

    Permutation& operator=(const Permutation& other) = default;
    Permutation& operator=(Permutation&& other) = default;
    Permutation& operator=(const std::initializer_list<size_t>& data) {
        size_t i = 0;
        for (auto it = data.begin(); it != data.end(); ++i, ++it) {
            data_[i] = *it;
        }
        return *this;
    }
    Permutation& operator=(const std::initializer_list<std::initializer_list<size_t>>& cycles) {
        std::iota(data_.begin(), data_.end(), 0);
        for (const auto& cycle : cycles) {
            size_t start = *cycle.begin();
            size_t prev = start;
            for (const auto& elem : cycle) {
                data_[prev] = elem;
                prev = elem;
            }
            data_[prev] = start;
        }
        return *this;
    }

    Permutation operator*(const Permutation& other) const {
        Permutation ans(size());
        for (size_t i = 0; i < size(); ++i) {
            ans.data_[i] = data_[other.data_[i]];
        }
        return ans;
    }
    Permutation& operator*=(const Permutation& other) {
        return *this = *this * other;
    }

    bool operator==(const Permutation& other) const = default;
    bool operator!=(const Permutation& other) const = default;

    Permutation Inverse() const {
        Permutation ans(size());
        for (size_t i = 0; i < size(); ++i) {
            ans.data_[data_[i]] = i;
        }
        return ans;
    }
    Permutation Power(size_t indicator) const {
        if (indicator == 0) return Permutation(size());
        if (indicator % 2 == 0) return (*this * *this).Power(indicator / 2);
        return Power(indicator - 1) * *this;
    }

    size_t operator[](size_t ind) const {
        return data_[ind];
    }

    size_t size() const {
        return data_.size();
    }

    int sign() const {
        std::vector<bool> used(size(), false);
        size_t cmp_comps = 0;
        for (size_t i = 0; i < size(); ++i) {
            if (!used[i]) {
                size_t vertex = i;
                while (!used[vertex]) {
                    used[vertex] = true;
                    vertex = data_[vertex];
                }
                ++cmp_comps;
            }
        }
        return ((size() - cmp_comps) % 2 == 0 ? 1 : -1);
    }

    friend class AllPermutations;

private:
    std::vector<size_t> data_;
};

template<size_t n>
std::ostream& operator<<(std::ostream& stream, const Permutation& p) {
    for (size_t i = 0; i < p.size(); ++i) {
        stream << p[i] + 1 << ' ';
    }
    return stream;
}

class AllPermutations {
public:
    class Iterator {
    public:
        Iterator(size_t n, bool is_end = false): p_(n), is_end_(is_end), n_(n) {}

        Iterator operator++() {
            for (ssize_t i = n_ - 2; i >= 0; --i) {
                if (p_.data_[i + 1] > p_.data_[i]) {
                    size_t to_swap = n_ - 1;
                    for (size_t j = i; j < n_; ++j) {
                        if (p_.data_[j] > p_.data_[i]) {
                            to_swap = j;
                        }
                    }
                    std::swap(p_.data_[to_swap], p_.data_[i]);
                    std::reverse(p_.data_.begin() + i + 1, p_.data_.end());
                    return *this;
                }
            }
            is_end_ = true;
            return *this;
        }
        Iterator operator++(int) {
            auto cp = *this;
            ++*this;
            return cp;
        }

        const Permutation& operator*() const {
            return p_;
        }

        bool operator==(const Iterator& other) const {
            if (other.is_end_ && is_end_) return true;
            if (other.is_end_ || is_end_) return false;
            return p_ == other.p_;
        }
        bool operator!=(const Iterator& other) const {
            return !operator==(other);
        }

    private:
        Permutation p_;
        bool is_end_;
        size_t n_;
    };

    AllPermutations(size_t n): n_(n) {}

    Iterator begin() const {
        return Iterator(n_, false);
    }
    Iterator end() const {
        return Iterator(n_, true);
    }
private:
    size_t n_;
};
