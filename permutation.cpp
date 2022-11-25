#include "mymath.h"
#include "permutation.h"

#include <algorithm>
#include <numeric>
#include <ostream>

Permutation::Permutation(size_t n): data(n) {
    std::iota(data.begin(), data.end(), 0);
}

Permutation::Permutation(const std::vector<size_t>& data): data(data) {}
Permutation::Permutation(const std::initializer_list<size_t>& data): data(data) {}
Permutation::Permutation(const std::initializer_list<std::initializer_list<size_t>>& cycles) {
    size_t n = 0;
    for (const auto& cycle : cycles) {
        n = std::max(n, *std::max_element(cycle.begin(), cycle.end()));
    }
    data.resize(n + 1);
    *this = cycles;
}

Permutation& Permutation::operator=(const std::initializer_list<size_t>& data) {
    this -> data = data;
    return *this;
}

Permutation& Permutation::operator=(const std::initializer_list<std::initializer_list<size_t>>& cycles) {
    std::iota(data.begin(), data.end(), 0);
    for (const auto& cycle : cycles) {
        size_t start = *cycle.begin();
        size_t prev = start;
        for (const auto& elem : cycle) {
            data[prev] = elem;
            prev = elem;
        }
        data[prev] = start;
    }
    return *this;
}

Permutation Permutation::operator*(const Permutation& other) const {
    Permutation ans(std::max(data.size(), other.data.size()));
    for (size_t i = 0; i < std::max(data.size(), other.data.size()); ++i) {
        size_t j = i;
        if (other.data.size() > j) {
            j = other.data[j];
        }
        if (data.size() > j) {
            j = data[j];
        }
        ans.data[i] = j;
    }
    return ans;
}

Permutation& Permutation::operator*=(const Permutation& other) {
    return *this = *this * other;
}

Permutation Permutation::Inverse() const {
    Permutation ans(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        ans.data[data[i]] = i;
    }
    return ans;
}

Permutation Permutation::Power(size_t indicator) const {
    return fastpow(*this, indicator);
}

size_t Permutation::operator[](size_t ind) const {
    return data[ind];
}

size_t Permutation::size() const {
    return data.size();
}

std::ostream& operator<<(std::ostream& stream, const Permutation& p) {
    for (size_t i = 0; i < p.size(); ++i) {
        stream << p[i] + 1 << ' ';
    }
    return stream;
}

AllPermutations::AllPermutations(size_t n): size(n) {}

AllPermutations::Iterator::Iterator(size_t size, bool is_end): p(size), is_end(is_end) {}

AllPermutations::Iterator AllPermutations::Iterator::operator++() {
    size_t n = p.size();
    for (ssize_t i = n - 2; i >= 0; --i) {
        if (p.data[i + 1] > p.data[i]) {
            size_t to_swap = n - 1;
            for (size_t j = i; j < n; ++j) {
                if (p.data[j] > p.data[i]) {
                    to_swap = j;
                }
            }
            std::swap(p.data[to_swap], p.data[i]);
            std::reverse(p.data.begin() + i + 1, p.data.end());
            return *this;
        }
    }
    is_end = true;
    return *this;
}

AllPermutations::Iterator AllPermutations::Iterator::operator++(int) {
    auto cp = *this;
    ++*this;
    return cp;
}

const Permutation& AllPermutations::Iterator::operator*() const {
    return p;
}

bool AllPermutations::Iterator::operator==(const Iterator& other) const {
    if (other.is_end && is_end) return true;
    if (other.is_end || is_end) return false;
    return p == other.p;
}

bool AllPermutations::Iterator::operator!=(const Iterator& other) const {
    return !operator==(other);
}

AllPermutations::Iterator AllPermutations::begin() const {
    return Iterator(size, false);
}

AllPermutations::Iterator AllPermutations::end() const {
    return Iterator(size, true);
}
