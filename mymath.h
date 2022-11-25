#pragma once

#include <cstddef>

template<typename T>
T myabs(const T& a) {
    return (a >= static_cast<T>(0) ? a : -a);
}

template<typename T>
T fastpow(const T& a, size_t n) {
    if (n == 0) return T(1);
    if (n % 2 == 0) {
        return fastpow(a * a, n / 2);
    }
    return fastpow(a, n - 1) * a;
}


