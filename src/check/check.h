#pragma once

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

template <typename T>
inline bool equal_vectors(const std::vector<T> &a, const std::vector<T> &b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] != b[i]) {
            return false;
        }
    }
    return true;
}

template <typename T>
inline bool equal_vectors_tol(const std::vector<T> &a, const std::vector<T> &b, T tolerance, bool verbose = false) {
    if (a.size() != b.size()) {
        if (verbose) {
            std::cout << "a.size() != b.size()" << std::endl;
        }
        return false;
    }
    for (size_t i = 0; i < a.size(); i++) {
        if (fabs(double(a[i]) - double(b[i])) > tolerance) {
            if (verbose) {
                double diff = fabs(double(a[i]) - double(b[i]));
                std::cout << a[i] << " != " << b[i] << " diff = " << diff << std::endl;
            }
            return false;
        }
    }
    return true;
}

template <typename T>
inline bool equal_scalars_tol(T a, T b, T tolerance) {
    if (fabs(double(a) - double(b)) > tolerance) {
        return false;
    }
    return true;
}

template <typename T>
inline bool equal_complex_vectors_tol(const std::vector<std::complex<T>> &a,
                                      const std::vector<std::complex<T>> &b,
                                      T tolerance_real,
                                      T tolerance_imag) {
    if (a.size() != b.size()) {
        return false;
    }
    for (size_t i = 0; i < a.size(); i++) {
        if (fabs(double(a[i].real()) - double(b[i].real())) > tolerance_real) {
            return false;
        }
        if (fabs(double(a[i].imag()) - double(b[i].imag())) > tolerance_imag) {
            return false;
        }
    }
    return true;
}
