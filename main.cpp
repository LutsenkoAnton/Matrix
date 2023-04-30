#include "fraction.h"
#include "float.h"
#include "integer.h"
#include "integer_mod.h"
#include "linear_operator.h"
#include "matrix.h"
#include "vector.h"
#include "vector_space.h"

#include <cassert>
#include <iostream>

using namespace std;

using Frac = Fraction<Integer>;
// using Frac = Float;

int solve1() {
    // Matrix<IntegerMod> a = {
    //     {-5_im, -3_im,  1_im, -3_im},
    //     { 5_im,  4_im, -1_im,  2_im},
    //     {-5_im, -2_im,  3_im, -2_im},
    //     { 9_im,  4_im, -1_im,  6_im}};
    // Matrix<IntegerMod> b = {
    //     {-1_im, -5_im,  1_im,  3_im},
    //     { 3_im,  8_im, -1_im, -4_im},
    //     {-3_im, -4_im,  3_im,  2_im},
    //     { 3_im,  6_im, -1_im, -2_im}};
    // Matrix<IntegerMod> c = {
    //     { 1_im, -4_im, -2_im, -3_im},
    //     {-1_im, -4_im, -3_im, -4_im},
    //     { 1_im,  8_im,  6_im,  5_im},
    //     { 1_im,  4_im,  2_im,  5_im}};
    // cerr << a << endl;
    Matrix<Frac> a = {
        {-5_fi, -3_fi,  1_fi, -3_fi},
        { 5_fi,  4_fi, -1_fi,  2_fi},
        {-5_fi, -2_fi,  3_fi, -2_fi},
        { 9_fi,  4_fi, -1_fi,  6_fi}};
    Matrix<Frac> b = {
        {-1_fi, -5_fi,  1_fi,  3_fi},
        { 3_fi,  8_fi, -1_fi, -4_fi},
        {-3_fi, -4_fi,  3_fi,  2_fi},
        { 3_fi,  6_fi, -1_fi, -2_fi}};
    Matrix<Frac> c = {
        { 1_fi, -4_fi, -2_fi, -3_fi},
        {-1_fi, -4_fi, -3_fi, -4_fi},
        { 1_fi,  8_fi,  6_fi,  5_fi},
        { 1_fi,  4_fi,  2_fi,  5_fi}};
    cout << LinearOperator<Frac>(a).GetJNF() << endl;
    cout << LinearOperator<Frac>(b).GetJNF() << endl;
    cout << LinearOperator<Frac>(c).GetJNF() << endl;
    /* We see that Jordan normal forms of A and B are equal, but C is not, so
     * A and B are the same operators
     * C is different from A and B
     */
    // cerr << Poly<Frac>({0_fi, 1_fi}) / Poly<Frac>::ONE() << endl;
    return 0;
}

int solve2() {
    Matrix<Frac> phi = {
        { 4_fi,  3_fi, -1_fi,  0_fi,  0_fi, 0_fi},
        { 0_fi,  6_fi, -1_fi,  0_fi,  0_fi, 0_fi},
        { 0_fi,  2_fi,  3_fi,  0_fi,  0_fi, 0_fi},
        {-3_fi, -2_fi,  4_fi,  2_fi, -1_fi, 0_fi},
        { 6_fi,  1_fi, -8_fi,  4_fi,  6_fi, 0_fi},
        {-3_fi, -1_fi,  4_fi, -2_fi, -2_fi, 5_fi}};
    auto J = LinearOperator<Frac>(phi).GetJNF();
    cout << "Jordan normal form:\n";
    cout << J << endl;
    auto C = LinearOperator<Frac>(phi).GetJordanBasis();
    C *= 18_fi; // Make numbers integer
    cout << "Jordan basis:\n";
    cout << C << endl;
    assert(C.Inverse() * phi * C == J);
    return 0;
}

int solve3() {
    Matrix<Frac> G = {
        {-1_fi,  2_fi, -1_fi},
        { 1_fi, -1_fi, -1_fi},
        { 2_fi, -3_fi,  1_fi}};
    Matrix<Frac> A = {
        {-4_fi,  8_fi, -8_fi},
        {-3_fi,  6_fi, -6_fi},
        { 1_fi, -1_fi,  0_fi}};
    Matrix<Frac> F = {
        {1_fi, -1_fi,  2_fi},
        {1_fi, -2_fi,  3_fi},
        {1_fi,  1_fi, -1_fi}};
    Matrix<Frac> B = {
        {-4_fi,  0_fi, -1_fi},
        { 8_fi, -4_fi,  8_fi},
        { 8_fi, -4_fi,  8_fi}};
    LinearOperator<Frac> phi(G * A * G.Inverse());
    LinearOperator<Frac> psi(F * B * F.Inverse());
    cout << phi.GetData() << endl;
    cout << psi.GetData() << endl;
    cout << "a)" << endl << (phi * psi).GetData() << endl;
    cout << "b)\n";
    auto sum = phi.ker() + psi.ker();
    for (const auto& vec : sum.GetBasis()) {
        cout << vec << endl;
    }
    cout << "\nc)\n";
    auto inter = phi.Im().inter(psi.Im());
    for (const auto& vec : inter.GetBasis()) {
        cout << vec << endl;
    }
    return 0;
}
int solve4() {
    Poly<Frac> t({0_fi, 1_fi});
    Matrix<Poly<Frac>> A = {
        {    2_fi,     1_fi,  0_fi,  0_fi,  0_fi},
        {    0_fi,     2_fi,  0_fi,  0_fi,  0_fi},
        {t + 1_fi,     0_fi, -2_fi,  0_fi,    t},
        {    0_fi, t - 1_fi,  1_fi, -2_fi,  0_fi},
        {    0_fi,     0_fi,  0_fi,  0_fi, -2_fi}};
    auto cp = A.CharPoly();
    for (auto [power, elem] : cp.GetCoefficients()) {
        assert(elem.deg() == 0);
    }
    // As we see all coefficients are integers so let's convert the polinomial to numbers
    Poly<Frac> numcp;
    for (auto [power, elem] : cp.GetCoefficients()) {
        numcp += Poly<Frac>(elem.GetCoefficients()[0], power);
    }
    std::vector<Frac> roots;
    for (auto [root, cnt] : numcp.try_solve()) {
        roots.push_back(root);
    }
    auto e = Matrix<Poly<Frac>>::IdentityMatrix(5);
    // Commented out because it produces too much output
    /* for (const auto& root : roots) {
        auto C = A - Poly<Frac>(root) * e;
        cout << C << endl;
        cout << C.Power(2) << endl;
        cout << C.Power(3) << endl;
        cout << C.Power(4) << endl;
    } */
    // From this we can see that rank of (A - 2E) is 4
    // (A - 2E)^2 is 3
    // (A - 2E)^2 is 3
    //
    // rank of (A + 2E) is different for x = 0 and x ≠ 0
    // same for (A + 2E)^2
    // ranks of (A + 2e)^3 and (A + 2E)^4 are 2
    // So, there are only two cases: x = 0 and x ≠ 0

    Matrix<Frac> B(5, 5);
    Matrix<Frac> C(5, 5);
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            B[i][j] = A[i][j](0_fi);
            C[i][j] = A[i][j](1_fi);
        }
    }
    cout << "x = 0:\n" << LinearOperator<Frac>(B).GetJNF() << endl;
    cout << "x ≠ 0:\n" << LinearOperator<Frac>(C).GetJNF() << endl;
    // See tex
    return 0;
}
int solve5() {
    // I am too lazy to write BilinearForm class sorry(
    Matrix<Frac> B = {
        { 0_fi,  0_fi,  2_fi, -3_fi},
        { 0_fi,  0_fi, -3_fi,  5_fi},
        { 2_fi, -3_fi,  2_fi, -4_fi},
        {-3_fi,  5_fi, -4_fi,  5_fi}};
    // Applying Jackobian method for reversed basis
    assert(Det(B) != 0_fi);
    assert(Det(B.Slice(3, 3, 1, 1)) != 0_fi);
    assert(Det(B.Slice(2, 2, 2, 2)) != 0_fi);
    assert(Det(B.Slice(1, 1, 3, 3)) != 0_fi); 
    Vector<Frac> v = {1_fi, 1_fi, 1_fi, 1_fi};
    VectorSpace<Frac> sp(std::vector({v}));
    sp.MakeFullBasis();
    auto basis = sp.GetBasis();
    reverse(basis.begin(), basis.end());
    std::vector<Vector<Frac>> new_basis(4, Vector<Frac>(4));
    for (int i = 3; i >= 0; --i) {
        new_basis[i] = basis[i];
        for (int j = i + 1; j < 4; ++j) {
            Frac beta1 = (static_cast<Matrix<Frac>>(basis[i]).Transpose() * B * new_basis[j])[0][0];
            Frac beta2 = (static_cast<Matrix<Frac>>(new_basis[j]).Transpose() * B * new_basis[j])[0][0];
            // cerr << basis[i] << ' ' << new_basis[j] << ' ' << beta1 << ' ' << beta2 << ' ' << i << ' '  << j << endl;
            new_basis[i] -= new_basis[j] * (beta1 / beta2);
        }
    }
    new_basis[0] *= 19_fi; // Making basis integer
    new_basis[1] *= 4_fi;
    auto C = VectorSpace(new_basis).GetBasisMatrix();
    cout << "New basis:\n";
    for (const auto& elem : new_basis) {
        cout << elem << endl;
    }
    cout << "Matrix in this basis:\n";
    cout << C.Transpose() * B * C;
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
