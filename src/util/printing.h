#pragma once

#include <iomanip>  // for setw
#include <iostream> // for cout, endl
#include <string>
#include <vector>

/// @brief Prints the components of a vector
/// @param T number type
/// @param name name of the vector, e.g., "x"
/// @param x vector data
template <typename T>
void print_vector(std::string name, const std::vector<T> &x) {
    std::cout << name << " = ";
    for (auto v : x) {
        std::cout << v << ", ";
    }
    std::cout << "\n";
}

/// @brief Prints a number in scientific format
/// @param T must be a number, e.g., float or double
/// @param value the number
/// @param width width of the whole output, e.g., 23
/// @param precision precision, e.g., 15
template <typename T>
void print_scientific(T value, size_t width, size_t precision) {
    std::cout << std::setw(width) << std::scientific << std::setprecision(precision) << value;
    std::cout << std::resetiosflags(std::ios_base::scientific) << std::resetiosflags(std::ios_base::fixed) << std::resetiosflags(std::ios_base::boolalpha);
}
