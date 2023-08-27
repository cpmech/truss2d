#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <complex>
#include <string>
#include <vector>

#include "../util/doctest.h"
#include "check.h"
using namespace std;

TEST_CASE("equal_vectors") {
    SUBCASE("int vectors") {
        vector<int> a{1, 2, 3};
        vector<int> b{1, 2, 3};
        vector<int> c{1, 2, 4};
        vector<int> d{1, 2};
        CHECK(equal_vectors(a, b) == true);
        CHECK(equal_vectors(b, c) == false);
        CHECK(equal_vectors(c, d) == false);
    }

    SUBCASE("size_t vectors") {
        vector<size_t> a{1, 2, 3};
        vector<size_t> b{1, 2, 3};
        vector<size_t> c{1, 2, 4};
        vector<size_t> d{1, 2};
        CHECK(equal_vectors(a, b) == true);
        CHECK(equal_vectors(b, c) == false);
        CHECK(equal_vectors(c, d) == false);
    }

    SUBCASE("string vectors") {
        vector<string> a{"1", "2", "3"};
        vector<string> b{"1", "2", "3"};
        vector<string> c{"1", "2", "4"};
        vector<string> d{"1", "2"};
        CHECK(equal_vectors(a, b) == true);
        CHECK(equal_vectors(b, c) == false);
        CHECK(equal_vectors(c, d) == false);
    }
}

TEST_CASE("equal_vectors_tol") {
    SUBCASE("float vectors") {
        float tolerance = 1e-6;
        vector<float> a{1.0f, 2.0f, 3.0f};
        vector<float> b{1.0f, 2.0f, 3.0f + tolerance};
        vector<float> c{1.0f, 2.0f, 3.0f + 1e-5};
        vector<float> d{1.0f, 2.0f};
        CHECK(equal_vectors_tol(a, b, tolerance) == true);
        CHECK(equal_vectors_tol(b, c, tolerance) == false);
        CHECK(equal_vectors_tol(c, d, tolerance) == false);
    }

    SUBCASE("double vectors") {
        double tolerance = 1e-15;
        vector<double> a{1.0, 2.0, 3.0};
        vector<double> b{1.0, 2.0, 3.0 + tolerance};
        vector<double> c{1.0, 2.0, 3.0 + 1e-14};
        vector<double> d{1.0, 2.0};
        CHECK(equal_vectors_tol(a, b, tolerance) == true);
        CHECK(equal_vectors_tol(b, c, tolerance) == false);
        CHECK(equal_vectors_tol(c, d, tolerance) == false);
    }
}

TEST_CASE("equal_complex_vectors_tol") {
    SUBCASE("float complex vectors") {
        float tol_real = 1e-6;
        float tol_imag = 1e-6;
        vector<complex<float>> a{1.0f, 2.0f, 3.0f};
        vector<complex<float>> b{1.0f, 2.0f, 3.0f + tol_real};
        vector<complex<float>> c{1.0f, 2.0f, 3.0f + 1e-5};
        vector<complex<float>> d{1.0f, 2.0f};
        vector<complex<float>> aa{complex<float>(1, 1), complex<float>(2, 2), complex<float>(3, 3)};
        vector<complex<float>> bb1{complex<float>(1, 1), complex<float>(2, 2), complex<float>(3 + tol_real, 3)};
        vector<complex<float>> bb2{complex<float>(1, 1), complex<float>(2, 2), complex<float>(3, 3 + tol_imag)};
        vector<complex<float>> bb3{complex<float>(1, 1), complex<float>(2, 2), complex<float>(3 + tol_real, 3 + tol_imag)};
        vector<complex<float>> cc{complex<float>(1, 1), complex<float>(2, 2), complex<float>(3, 3 + 1e-5)};
        vector<complex<float>> dd{complex<float>(1, 1), complex<float>(2, 2)};
        CHECK(equal_complex_vectors_tol(a, b, tol_real, tol_imag) == true);
        CHECK(equal_complex_vectors_tol(b, c, tol_real, tol_imag) == false);
        CHECK(equal_complex_vectors_tol(c, d, tol_real, tol_imag) == false);
        CHECK(equal_complex_vectors_tol(aa, bb1, tol_real, tol_imag) == true);
        CHECK(equal_complex_vectors_tol(aa, bb2, tol_real, tol_imag) == true);
        CHECK(equal_complex_vectors_tol(aa, bb3, tol_real, tol_imag) == true);
        CHECK(equal_complex_vectors_tol(bb1, cc, tol_real, tol_imag) == false);
        CHECK(equal_complex_vectors_tol(bb2, cc, tol_real, tol_imag) == false);
        CHECK(equal_complex_vectors_tol(bb3, cc, tol_real, tol_imag) == false);
        CHECK(equal_complex_vectors_tol(cc, dd, tol_real, tol_imag) == false);
    }

    SUBCASE("double complex vectors") {
        double tol_real = 1e-15;
        double tol_imag = 1e-15;
        vector<complex<double>> a{1.0, 2.0, 3.0};
        vector<complex<double>> b{1.0, 2.0, 3.0 + tol_real};
        vector<complex<double>> c{1.0, 2.0, 3.0 + 1e-14};
        vector<complex<double>> d{1.0, 2.0};
        vector<complex<double>> aa{complex<double>(1, 1), complex<double>(2, 2), complex<double>(3, 3)};
        vector<complex<double>> bb1{complex<double>(1, 1), complex<double>(2, 2), complex<double>(3 + tol_real, 3)};
        vector<complex<double>> bb2{complex<double>(1, 1), complex<double>(2, 2), complex<double>(3, 3 + tol_imag)};
        vector<complex<double>> bb3{complex<double>(1, 1), complex<double>(2, 2), complex<double>(3 + tol_real, 3 + tol_imag)};
        vector<complex<double>> cc{complex<double>(1, 1), complex<double>(2, 2), complex<double>(3, 3 + 1e-14)};
        vector<complex<double>> dd{complex<double>(1, 1), complex<double>(2, 2)};
        CHECK(equal_complex_vectors_tol(a, b, tol_real, tol_imag) == true);
        CHECK(equal_complex_vectors_tol(b, c, tol_real, tol_imag) == false);
        CHECK(equal_complex_vectors_tol(c, d, tol_real, tol_imag) == false);
        CHECK(equal_complex_vectors_tol(aa, bb1, tol_real, tol_imag) == true);
        CHECK(equal_complex_vectors_tol(aa, bb2, tol_real, tol_imag) == true);
        CHECK(equal_complex_vectors_tol(aa, bb3, tol_real, tol_imag) == true);
        CHECK(equal_complex_vectors_tol(bb1, cc, tol_real, tol_imag) == false);
        CHECK(equal_complex_vectors_tol(bb2, cc, tol_real, tol_imag) == false);
        CHECK(equal_complex_vectors_tol(bb3, cc, tol_real, tol_imag) == false);
        CHECK(equal_complex_vectors_tol(cc, dd, tol_real, tol_imag) == false);
    }
}
