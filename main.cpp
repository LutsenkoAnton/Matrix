#include <complex>
#include <iostream>
#include <iomanip>
#include <type_traits>

#include "fraction.h"
#include "matrix.h"
#include "permutation.h"
#include "rational_function.h"

using namespace std;

int solve1() {
    // SquareMatrix<int> a = 
       // {{0, 0, 0, 1, 0, 0, 0, 0},
        // {0, 0, 0, 0, 0, 0, 1, 0},
        // {0, 0, 1, 0, 0, 0, 0, 0},
        // {0, 0, 0, 0, 1, 0, 0, 0},
        // {0, 0, 0, 0, 0, 0, 0, 1},
        // {1, 0, 0, 0, 0, 0, 0, 0},
        // {0, 0, 0, 0, 0, 1, 0, 0},
        // {0, 1, 0, 0, 0, 0, 0, 0}};
    // SquareMatrix<int> b = 
       // {{0, 0, 0, 0, 0, 1, 0, 0},
        // {0, 0, 0, 0, 0, 0, 1, 0},
        // {0, 1, 0, 0, 0, 0, 0, 0},
        // {0, 0, 0, 0, 1, 0, 0, 0},
        // {0, 0, 1, 0, 0, 0, 0, 0},
        // {0, 0, 0, 1, 0, 0, 0, 0},
        // {1, 0, 0, 0, 0, 0, 0, 0},
        // {0, 0, 0, 0, 0, 0, 0, 1}};
    // SquareMatrix<int> c = 
       // {{0, 0, 0, 0, 0, 0, 1, 0},
        // {0, 0, 0, 0, 0, 1, 0, 0},
        // {0, 1, 0, 0, 0, 0, 0, 0},
        // {0, 0, 0, 1, 0, 0, 0, 0},
        // {0, 0, 1, 0, 0, 0, 0, 0},
        // {0, 0, 0, 0, 0, 0, 0, 1},
        // {1, 0, 0, 0, 0, 0, 0, 0},
        // {0, 0, 0, 0, 1, 0, 0, 0}};
    Permutation a = {3, 6, 2, 4, 7, 0, 5, 1};
    Permutation b = {5, 6, 1, 4, 2, 3, 0, 7};
    Permutation c = {{0, 6}, {1, 5, 7, 4, 2, 3}};
    // s * c * s = (a ^ -1 * b ^ 17) ^ 161
    auto x = (a.Inverse() * b.Power(17)).Power(161);
    for(const auto& s : AllPermutations(8)) {
        if (s * c * s == x) {
            cout << s << endl;
        }
    }
    return 0;
}

int solve2() {
    SquareMatrix<Fraction> a = 
       {{3, 1, -2, 1},
        {-1, -2, 1, 3},
        {-3, -2, 1, 1},
        {3, -1, -1, -1}};
    SquareMatrix<Fraction> b = 
       {{-3, -3, -3, -3},
        {2, 1, -2, -3},
        {1, -1, 2, 1},
        {2, -1, -2, -3}};
    SquareMatrix<Fraction> c = 
       {{1, 1, 1, -1},
        {-1, -2, -2, 2},
        {1, 2, 3, -3},
        {-1, -2, -3, 2}};
    SquareMatrix<Fraction> d = 
       {{10, -8, -3, -2},
        {10, -6, 7, -3},
        {-8, 3, 3, 6},
        {5, -4, 9, -4}};
    // a * (x + d) ^ -1 * b = c
    // (x + d)^-1 = a^-1 * c * b^-1
    // x = (a ^ -1 * c * b^-1) - d;
    auto x = (a.Inverse() * c * b.Inverse()).Inverse() - d;
    assert(a * (x + d).Inverse() * b == c);
    cout << x << endl;
    return 0;
}

int solve3() {
    SquareMatrix<Fraction> a = 
       {{-1, 4, -2, -4},
        {1, -2, -2, -5},
        {-1, 2, 2, 5},
        {0, 0, -2, -5}};
    auto p = a.CharPoly();
    cout << p << endl;
    // det((a^2 - 1)^-2) = det((a - 1)(a + 1)) ^ -2 = (det(a - 1) * det(a + 1))^-2 = (det(1 - a) * det(-1 - a))^-2 = 1 / (p(1) ^ 2 * p(-1) ^ 2)
    cout << Fraction(1)/(p(1) * p(1) * p(-1) * p(-1));
    cout << endl;
    return 0;
}

int solve4() {
    Poly x({0LL, 1LL});
    SquareMatrix<RationalFunction<int64_t>> a = 
       {{x, 4, -7, -7, -9, -8, x},
        {-4, x, -2, 5, -6, 2, 4},
        {-1, -5, x, 1, 5, -5, 2},
        {-4, x, 5, 1, -3, -9, -6},
        {8, 6, -7, -2, x, 7, -5},
        {1, -9, -4, -3, 2, x, -1},
        {-4, -6, -2, x, 2, 7, -1}};
    cout << PermDet(a).force_devide() << endl;
    // cout << Det(a).force_devide() << endl; // Works faster, but overflows
    return 0;
}

int solve5() {
    Matrix<int64_t> a =
       {{1, -1},
        {-1, 1},
        {4, -4},
        {-3, 3},
        {2, -2}};
    Matrix<int64_t> b =
       {{-3, -1, 2, 1, 3},
        {2, 3, -2, -3, 1}};
    cout << static_cast<SquareMatrix<int64_t>>(a * b).CharPoly() << endl;
    cout << static_cast<SquareMatrix<int64_t>>(a * b).PermCharPoly() << endl;
    return 0;
}

int main() {
    // I don't regret a single minute spent on this thing
    solve1();
    cout << "--------------------------------------\n";
    solve2();
    cout << "--------------------------------------\n";
    solve3();
    cout << "--------------------------------------\n";
    solve4();
    cout << "--------------------------------------\n";
    solve5();
    return 0;
}
