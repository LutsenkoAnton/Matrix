#pragma once

#include "variable.h"

#include <iostream>
#include <map>

template<typename T>
class LinearCombination {
public:
    LinearCombination() = default;
    LinearCombination(const LinearCombination& other) = default;
    LinearCombination(LinearCombination&& other) = default;
    LinearCombination& operator=(const LinearCombination& other) = default;
    LinearCombination& operator=(LinearCombination&& other) = default;

    LinearCombination(Variable v): coefficients_({{v, 1}}) {}
    LinearCombination(T coefficient): coefficients_({{Variable(0), coefficient}}) {}
    LinearCombination(std::map<Variable, T> coefficients): coefficients_(coefficients) {}

    LinearCombination operator+(const LinearCombination& other) const {
        LinearCombination ans;
        for (const auto& [var, coefficient] : coefficients_) {
            ans.coefficients_[var] += coefficient;
        }
        for (const auto& [var, coefficient] : other.coefficients_) {
            ans.coefficients_[var] += coefficient;
        }
        ans.Clear();
        return ans;
    }
    LinearCombination operator-() const {
        LinearCombination ans(coefficients_);
        for (auto &[var, coefficient]: ans) {
            coefficient = -coefficient;
        }
        ans.Clear();
        return ans;
    }
    LinearCombination operator-(const LinearCombination& other) const {
        return *this + -other;
    }
    LinearCombination operator*(const T& other) const {
        LinearCombination ans(coefficients_);
        for (auto &[var, coefficient]: ans) {
            coefficient = coefficient * other;
        }
        ans.Clear();
        return ans;
    }
    LinearCombination& operator+=(const LinearCombination& other) {
        *this = *this + other;
        return *this;
    }
    LinearCombination& operator-=(const LinearCombination& other) {
        *this = *this - other;
        return *this;
    }
    LinearCombination& operator*=(const T& other) {
        *this = *this * other;
        return *this;
    }
    T& operator[](const Variable& var) {
        return coefficients_[var];
    }
    T operator[](const Variable& var) const {
        auto it = coefficients_.find(var);
        if (it == coefficients_.end()) return T(0);
        return it -> second;
    }

    bool is_constant() const {
        if (coefficients_.size() > 1) return false;
        if (coefficients_.size() == 0) return true;
        if (coefficients_.count(Variable(0)) == 1) return true;
        return false;
    }

    typename std::map<Variable, T>::iterator begin() {
        return coefficients_.begin();
    }

    typename std::map<Variable, T>::const_iterator begin() const {
        return coefficients_.begin();
    }

    typename std::map<Variable, T>::iterator end() {
        return coefficients_.end();
    }

    typename std::map<Variable, T>::const_iterator end() const {
        return coefficients_.end();
    }

private:
    void Clear() {
        for (auto it = coefficients_.begin(); it != coefficients_.end();) {
            if (it -> second == 0) {
                it = coefficients_.erase(it);
            } else {
                ++it;
            }
        }
    }
    std::map<Variable, T> coefficients_;
};

template<typename T>
LinearCombination<T> operator*(T coeff, LinearCombination<T> comb) {
    return comb * coeff;
}

template<typename T>
std::ostream& operator<<(std::ostream& stream, LinearCombination<T> comb) {
    T free_monom = 0;
    bool first = true;
    for (const auto &[var, coeff] : comb) {
        if (var.id == 0) {
            free_monom = coeff;
        } else {
            if (first) {
                first = false;
            } else {
                stream << " + ";
            }
            if (coeff == 1) {
                stream << "x[" << var.id << "]";
            } else {
                stream << coeff << " * x[" << var.id << "]";
            }
        }
    }
    if (free_monom != T(0)) {
        if (first) {
            stream << free_monom;
        } else {
            stream << " + " << free_monom;
        }
    }
    return stream;
}

