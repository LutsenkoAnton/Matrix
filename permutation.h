#pragma once

#include <array>
#include <initializer_list>
#include <iostream>
#include <numeric>

template<size_t n>
class Permutation {
public:
    Permutation() {
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
        Permutation ans;
        for (size_t i = 0; i < n; ++i) {
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
        Permutation ans;
        for (size_t i = 0; i < n; ++i) {
            ans.data_[data_[i]] = i;
        }
        return ans;
    }
    Permutation Power(size_t indicator) const {
        if (indicator == 0) return Permutation();
        if (indicator % 2 == 0) return (*this * *this).Power(indicator / 2);
        return Power(indicator - 1) * *this;
    }

    size_t operator[](size_t ind) const {
        return data_[ind];
    }

    size_t size() const {
        return n;
    }
    int sign() const {
        std::array<bool, n> used = {};
        size_t cmp_comps = 0;
        for (size_t i = 0; i < n; ++i) {
            if (!used[i]) {
                size_t vertex = i;
                while (!used[vertex]) {
                    used[vertex] = true;
                    vertex = data_[vertex];
                }
                ++cmp_comps;
            }
        }
        return ((n - cmp_comps) % 2 == 0 ? 1 : -1);
    }

    template<size_t n_>
    friend class AllPermutations;

    static inline const Permutation ZERO = Permutation(); // This is to satisfy RingWithOne
    static inline const Permutation ONE = Permutation();

private:
    std::array<size_t, n> data_;
};

template<size_t n>
std::ostream& operator<<(std::ostream& stream, const Permutation<n>& p) {
    for (size_t i = 0; i < p.size(); ++i) {
        stream << p[i] + 1 << ' ';
    }
    return stream;
}

template<size_t n>
class AllPermutations {
public:
    class Iterator {
    public:
        Iterator(bool is_end = false): is_end_(is_end) {}

        Iterator operator++() {
            for (ssize_t i = n - 2; i >= 0; --i) {
                if (p_.data_[i + 1] > p_.data_[i]) {
                    size_t to_swap = n - 1;
                    for (size_t j = i; j < n; ++j) {
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

        const Permutation<n>& operator*() const {
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
        Permutation<n> p_;
        bool is_end_;
    };

    Iterator begin() const {
        return Iterator(false);
    }
    Iterator end() const {
        return Iterator(true);
    }
};
