#pragma once

#include <exception>

class WrongSizeException : public std::exception {
public:
    const char* what() const noexcept override {
        return "Wrong sizes of matrixes";
    }
};

class SingularMatrixException : public std::exception {
public:
    const char* what() const noexcept override {
        return "Matrix is singular";
    }
};

