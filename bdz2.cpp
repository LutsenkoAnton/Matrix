#include <iostream>

#include "fraction.h"
#include "matrix.h"
#include "permutation.h"
#include "integer.h"

using Frac = Fraction<Integer>;
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
    Matrix<Frac, 4> a =
       {{3_fi, 1_fi, -2_fi, 1_fi},
        {-1_fi, -2_fi, 1_fi, 3_fi},
        {-3_fi, -2_fi, 1_fi, 1_fi},
        {3_fi, -1_fi, -1_fi, -1_fi}};
    Matrix<Frac, 4> b =
       {{-3_fi, -3_fi, -3_fi, -3_fi},
        {2_fi, 1_fi, -2_fi, -3_fi},
        {1_fi, -1_fi, 2_fi, 1_fi},
        {2_fi, -1_fi, -2_fi, -3_fi}};
    Matrix<Frac, 4> c =
       {{1_fi, 1_fi, 1_fi, -1_fi},
        {-1_fi, -2_fi, -2_fi, 2_fi},
        {1_fi, 2_fi, 3_fi, -3_fi},
        {-1_fi, -2_fi, -3_fi, 2_fi}};
    Matrix<Frac, 4> d =
       {{10_fi, -8_fi, -3_fi, -2_fi},
        {10_fi, -6_fi, 7_fi, -3_fi},
        {-8_fi, 3_fi, 3_fi, 6_fi},
        {5_fi, -4_fi, 9_fi, -4_fi}};
    // SquareMatrix<Frac> a = matrix_cast<Frac>(ia);
    // SquareMatrix<Frac> b = matrix_cast<Frac>(ib);
    // SquareMatrix<Frac> c = matrix_cast<Frac>(ic);
    // SquareMatrix<Frac> d = matrix_cast<Frac>(id);
    // a * (x + d) ^ -1 * b = c
    // (x + d)^-1 = a^-1 * c * b^-1
    // x = (a ^ -1 * c * b^-1)^-1 - d;
    // Throws an error if a or b are not invertible
    auto x = (a.Inverse() * c * b.Inverse()).Inverse() - d;
    assert(a * (x + d).Inverse() * b == c);
    cout << x << endl;
    return 0;
}

int solve3() {
    Matrix<Frac, 4> a =
       {{-1_fi, 4_fi, -2_fi, -4_fi},
        {1_fi, -2_fi, -2_fi, -5_fi},
        {-1_fi, 2_fi, 2_fi, 5_fi},
        {0_fi, 0_fi, -2_fi, -5_fi}};
    auto p = a.CharPoly();
    cout << p << endl;
    Integer o = 1;
    // det((a^2 - 1)^-2) = det((a - 1)(a + 1)) ^ -2 = (det(a - 1) * det(a + 1))^-2 = (det(1 - a) * det(-1 - a))^-2 = 1 / (p(1) ^ 2 * p(-1) ^ 2)
    cout << Frac(1)/(p(o) * p(o) * p(-o) * p(-o)) << endl;
    // cout << Det((a * a - SquareMatrix<Frac>::IdentityMatrix(4)).Inverse().Power(2));
    cout << endl;
    return 0;
}

int solve4() {
    Poly<Frac> x({0_fi, 1_fi});
    Matrix<Poly<Frac>, 7> a =
       {{x, 4_fi, -7_fi, -7_fi, -9_fi, -8_fi, x},
        {-4_fi, x, -2_fi, 5_fi, -6_fi, 2_fi, 4_fi},
        {-1_fi, -5_fi, x, 1_fi, 5_fi, -5_fi, 2_fi},
        {-4_fi, x, 5_fi, 1_fi, -3_fi, -9_fi, -6_fi},
        {8_fi, 6_fi, -7_fi, -2_fi, x, 7_fi, -5_fi},
        {1_fi, -9_fi, -4_fi, -3_fi, 2_fi, x, -1_fi},
        {-4_fi, -6_fi, -2_fi, x, 2_fi, 7_fi, -1_fi}};
    auto p = Det(a); // Owerflows, but gives correct result
    // auto p = PermDet(a).force_devide(); // Works slower
    cout << p << endl;
    cout << p.GetCoefficients()[5] << endl;
    return 0;
}

int solve5() {
    Matrix<Frac, 5, 2> a =
       {{1_fi, -1_fi},
        {-1_fi, 1_fi},
        {4_fi, -4_fi},
        {-3_fi, 3_fi},
        {2_fi, -2_fi}};
    Matrix<Frac, 2, 5> b =
       {{-3_fi, -1_fi, 2_fi, 1_fi, 3_fi},
        {2_fi, 3_fi, -2_fi, -3_fi, 1_fi}};
    cout << (a * b).CharPoly() << endl;
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
