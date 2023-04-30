#pragma once

//These concepts are made so that i can better understand what methods do i need

template<typename T>
concept MultiplicationMonoid = requires (T a, T b) {
    T::ONE();
    a * b;
};

template<typename T>
concept RingWithOne = requires (T a, T b) { // I can't check basic acsioms, but it is convenient to require these operators
    a + b;
    a - b;
    T::ZERO();
    T::ONE();
    -a;
    a * b;
};

template<typename T>
concept EuclideanRing = RingWithOne<T> && requires(T a, T b){
    a % b;
    a / b;
};

template<typename T>
concept Field = RingWithOne<T> && !EuclideanRing<T> && requires(T a, T b) {
    a / b;
};
