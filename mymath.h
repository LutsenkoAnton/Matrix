#pragma once

#include "myconcepts.h"
#include "poly.h"

#include <cstddef>
#include <iostream>

template<RingWithOne T>
T abs(const T& a) {
    return (a >= T::ZERO() ? a : -a);
}

template<MultiplicationMonoid T>
T fastpow(const T& a, size_t n) {
    if (n == 0) return T::ONE();
    if (n % 2 == 0) {
        return fastpow(a * a, n / 2);
    }
    return fastpow(a, n - 1) * a;
}

template<EuclideanRing T>
T gcd(const T& a, const T& b) {
    if (b == T::ZERO()) {
        return a;
    }
    return gcd(b, a % b);
}
