#include "fraction.h"
#include "integer.h"
#include "linear_operator.h"
#include "matrix.h"
#include "vector.h"
#include "vector_space.h"

#include <cassert>
#include <iostream>

using namespace std;

using Frac = Fraction<Integer>;

int solve1() {
    Matrix<Frac> a =
        {{-1_fi, -1_fi,  1_fi, -3_fi, -2_fi,  3_fi},
         {-1_fi, -1_fi,  2_fi, -5_fi, -1_fi,  4_fi},
         { 1_fi,  1_fi, -1_fi,  3_fi,  2_fi, -3_fi}};
    Vector<Frac> v1 = {-5_fi, 1_fi, -7_fi, -1_fi, 3_fi, 2_fi};
    Vector<Frac> v2 = {-8_fi, -1_fi, 0_fi, 2_fi, 3_fi, 1_fi};
    cerr << static_cast<Matrix<Frac>>(v1) << endl;
    assert((a * static_cast<Matrix<Frac>>(v1) == Matrix<Frac>(3, 1))); // Checking that Av1 = 0
    assert((a * static_cast<Matrix<Frac>>(v2) == Matrix<Frac>(3, 1))); // Checking that Av2 = 0
    auto basis = VectorSpace(a).GetBasis();
    std::vector<Vector<Frac>> new_basis = {v1, v2};
    for (const auto& vec : basis) {
        new_basis.push_back(vec);
    }
    // cout << new_basis.size() << endl;
    std::vector<Vector<Frac>> arr(6, Vector<Frac>(6)); // new_basis.size() == 4; +2 new vectors --- v1 and v2 // Yes, I have changed vectors to arrays because I forgot that sizes are known only in runtime
    for (size_t i = 0; i < 6; ++i) {
        arr[i] = new_basis[i];
    }
    auto ans = VectorSpace(arr).GetBasis(); // Basis will contain v1 and v2, beacause they are at the beginning of the array, so they can only be excluded by themselves
    for (const auto& vec : ans) {
        cout << vec << endl;
    }
    return 0;
}

int solve2() {
    Matrix<Frac> a = 
       {{-5_fi, -5_fi, -6_fi,  9_fi},
        { 5_fi,  5_fi,  6_fi, -9_fi},
        { 4_fi,  4_fi,  3_fi, -6_fi},
        {-4_fi, -4_fi, -3_fi,  6_fi}};
    cout << "Basis intersection:\n";
    LinearOperator phi(a);
    auto intersection = phi.Im().inter(phi.ker());
    for (const auto& vec : intersection.GetBasis()){
        cout << vec << endl;
    }
    cout << "Basis sum:\n";
    auto sum = phi.Im() + phi.ker();
    for (const auto& vec : sum.GetBasis()) {
        cout << vec << endl;
    }
    return 0;
}

int solve3() {
    Matrix<Frac> a =
       {{-2_fi, -2_fi,   7_fi, -4_fi},
        { 4_fi,  4_fi, -10_fi,  5_fi},
        { 0_fi,  0_fi,  -2_fi,  2_fi},
        { 0_fi,  0_fi,  -4_fi,  4_fi}};
    Vector<Frac> v = {9_fi, 1_fi, -1_fi, 5_fi};
    LinearOperator phi(a);
    LinearOperator psi(a.Transpose());
    VectorSpace im = phi.Power(4).Im(); // By Stabilization Lemma Power can be only 4
    VectorSpace ker = psi.Power(4).ker(); // By Stabilization Lemma Power can be only 4  
    assert((im.inter(ker).dim() == 0  && (im + ker).dim() == 4));
    Matrix<Frac> D(4, 2 * im.dim());
    Matrix<Frac> basisKer(4, 2 * im.dim());
    for (size_t i = 0; i < im.dim(); ++i) {
        for (size_t j = 0; j < 4; ++j) {
            D[j][i] = im.GetBasis()[i][j];
        }
    }
    for (size_t i = 0; i < im.dim(); ++i) {
        for (size_t j = 0; j < 4; ++j) {
            D[j][i + im.dim()] = ker.GetBasis()[i][j];
            basisKer[j][i + im.dim()] = ker.GetBasis()[i][j];
        }
    }
    Matrix<Frac> x = D.Inverse() * static_cast<Matrix<Frac>>(v);
    cout << "im part: " << Vector((D - basisKer) * x) << endl;
    cout << "ker part: " << Vector(basisKer * x) << endl;
    return 0;
}

int solve4() {
    Matrix<Integer> a = 
       {{-1_i, -2_i,  7_i,  4_i},
        { 4_i,  4_i, -6_i, -3_i},
        {-3_i, -2_i, -3_i, -2_i},
        { 6_i,  4_i,  6_i,  4_i}};

    cout << a.CharPoly() << endl; // x^4 - 4x^3 + 4x^2 = (x - 2) ^ 2 * x ^ 2
    Matrix<Frac> fa = 
       {{-1_fi, -2_fi,  7_fi,  4_fi},
        { 4_fi,  4_fi, -6_fi, -3_fi},
        {-3_fi, -2_fi, -3_fi, -2_fi},
        { 6_fi,  4_fi,  6_fi,  4_fi}};
    LinearOperator phi(fa);
    Frac l1 = 2_fi;
    Frac l2 = 0_fi;
    auto id = LinearOperator<Frac>::ONE(4);
    auto im = (phi - id * l1).Power(4).ker(); 
    im.GetBasis()[1] *= 2_fi; // For integer results; can be done without this line, but with it this is more beautiful
    auto ker = (phi - id * l2).Power(4).ker(); 
    LinearOperator psi(im, ker);
    cout << psi.GetData() << endl;
    return 0;
}

int solve5() {
    Matrix<Frac> G =
       {{-1_fi, -1_fi, 2_fi},
        {-1_fi, -2_fi, 1_fi},
        { 1_fi,  2_fi, 0_fi}};
    Matrix<Frac> v = {{-6_fi}, {-8_fi}, {6_fi}};
    Matrix<Frac> u = {{-8_fi}, {2_fi}, {-8_fi}};
    // C = G^-1 * A * G
    // A = G * C * G^-1
    // Av = u
    // G * C * G^-1 v = u
    // CG^-1v = G^-1u
    auto gv = G.Inverse() * v;
    auto gu = G.Inverse() * u;
    Matrix<Frac> C =
      {{ 0_fi,  1_fi,  6_fi},
       { 1_fi,  0_fi, -3_fi},
       {-2_fi, -1_fi,  0_fi}};
    C[0][0] = (gu[0][0] - C[0][1] * gv[1][0] - C[0][2] * gv[2][0]) / gv[0][0];
    C[1][1] = (gu[1][0] - C[1][0] * gv[0][0] - C[1][2] * gv[2][0]) / gv[1][0];
    C[2][2] = (gu[2][0] - C[2][0] * gv[0][0] - C[2][1] * gv[1][0]) / gv[2][0];
    auto A = G * C * G.Inverse();
    assert(A * v == u);
    cout << A;
    return 0;
}

int main() {
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
