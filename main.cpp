#include <iostream>
#include "matrix.h"

int main() {
    Matrix<int> a({{1, -4, -3}, {2, 1, 3}, {4, 3, 7}, {-4, 2, -2}, {-3, -2, -5}});
    Matrix<int> b({{1, 0, 0, 1}, {3, 0, 1, 0}, {-2, 1, 0, 0}});
    std::cout << a * b << std::endl;
}
