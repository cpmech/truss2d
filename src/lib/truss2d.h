#pragma once

#include "laclib.h"
#include <memory>
#include <vector>

/// @brief Implements a finite element solver for trusses in 2D
struct Truss2D {
    /// @brief Holds the number of nodes = coordinates.size() / 2
    size_t number_of_nodes;

    /// @brief Holds the number of elements = connectivity.size() / 2
    size_t number_of_elements;

    /// @brief Holds the total number of DOFs = 2 * number_of_nodes
    size_t total_ndof;

    /// @brief Coordinates x0 y0  x1 y1  ...  xnn ynn (size = 2 * number_of_nodes)
    std::vector<double> coordinates;

    /// @brief Connectivity 0 1  0 2  1 2  (size = 2 * number_of_elements)
    std::vector<size_t> connectivity;

    /// @brief Properties = E*A (size = number_of_elements)
    std::vector<double> properties;

    /// @brief Essential (displacement) prescribed?
    std::vector<bool> essential_prescribed;

    /// @brief Essential (displacement) boundary conditions (size = total_ndof)
    std::vector<double> essential_boundary_conditions;

    /// @brief Natural (force) boundary conditions (size = total_ndof)
    std::vector<double> natural_boundary_conditions;

    /// @brief Divide E*A by each length?
    bool divide_property_by_length;

    /// @brief Element stiffness matrix (4 x 4)
    std::unique_ptr<Matrix> kk_element;

    /// @brief Global displacements (size = total_ndof)
    std::vector<double> uu;

    /// @brief Global forces (size = total_ndof)
    std::vector<double> ff;

    /// @brief Internal forces of an element (size = total_ndof)
    std::vector<double> ff_int;

    /// @brief Residual = F - F_int (size = total_ndof)
    std::vector<double> residual;

    /// @brief Global displacements increments (size = total_ndof)
    std::vector<double> delta_uu;

    /// @brief Global forces increments (size = total_ndof)
    std::vector<double> delta_ff;

    /// @brief Internal forces increments of an element (size = total_ndof)
    std::vector<double> delta_ff_int;

    /// @brief Global stiffness matrix (nnz ~ 4 * 4 * number_of_elements)
    std::unique_ptr<CooMatrix> kk;

    /// @brief Allocates a new Truss2D structure
    /// @param coordinates x0 y0  x1 y1  ...  xnn ynn (size = 2 * number_of_nodes)
    /// @param connectivity 0 1  0 2  1 2  (size = 2 * number_of_elements)
    /// @param properties E*A (size = number_of_elements)
    /// @param essential_prescribed prescribed displacement flags
    /// @param essential_boundary_conditions prescribed displacement values
    /// @param natural_boundary_conditions prescribed external force values
    /// @param divide_property_by_length divide E*A by each rod length
    inline static std::unique_ptr<Truss2D> make_new(
        std::vector<double> coordinates,
        std::vector<size_t> connectivity,
        std::vector<double> properties,
        std::vector<bool> essential_prescribed,
        std::vector<double> essential_boundary_conditions,
        std::vector<double> natural_boundary_conditions,
        bool divide_property_by_length = true) {

        auto number_of_nodes = coordinates.size() / 2;
        auto number_of_elements = connectivity.size() / 2;
        auto total_ndof = 2 * number_of_nodes;
        auto nnz_max = 16 * number_of_elements; // can be optimized

        // StoredLayout layout = UPPER_TRIANGULAR;
        StoredLayout layout = FULL_MATRIX;

        return std::unique_ptr<Truss2D>{new Truss2D{
            number_of_nodes,
            number_of_elements,
            total_ndof,
            coordinates,
            connectivity,
            properties,
            essential_prescribed,
            essential_boundary_conditions,
            natural_boundary_conditions,
            divide_property_by_length,
            Matrix::make_new(4, 4),                           // kk_element
            std::vector<double>(total_ndof),                  // uu
            std::vector<double>(total_ndof),                  // ff
            std::vector<double>(total_ndof),                  // ff_int
            std::vector<double>(total_ndof),                  // residual
            std::vector<double>(total_ndof),                  // delta_uu
            std::vector<double>(total_ndof),                  // delta_ff
            std::vector<double>(total_ndof),                  // delta_ff_int
            CooMatrix::make_new(layout, total_ndof, nnz_max), // kk
        }};
    }

    /// @brief Returns the abscissa of node n (zero based)
    /// @param i index of node in 0 <= n < number_of_nodes
    inline double get_x(size_t n) const {
        if (n >= coordinates.size()) {
            throw "cannot get x coordinate because the node index is out-of-range";
        }
        return coordinates[n * 2];
    }

    /// @brief Returns the ordinate of node n (zero based)
    /// @param i index of node in 0 <= n < number_of_nodes
    inline double get_y(size_t n) const {
        if (n >= coordinates.size()) {
            throw "cannot get y coordinate because the node index is out-of-range";
        }
        return coordinates[n * 2 + 1];
    }

    /// @brief Returns the distance between node a and b
    /// @param i index of node in 0 <= a < number_of_nodes
    /// @param j index of node in 0 <= b < number_of_nodes
    inline double calculate_length(size_t a, size_t b) const {
        if (a >= coordinates.size()) {
            throw "cannot calculate the length because the a-th node index is out-of-range";
        }
        if (b >= coordinates.size()) {
            throw "cannot calculate the length because the b-th node index is out-of-range";
        }
        auto xa = coordinates[a * 2];
        auto ya = coordinates[a * 2 + 1];
        auto xb = coordinates[b * 2];
        auto yb = coordinates[b * 2 + 1];
        auto dx = xb - xa;
        auto dy = yb - ya;
        return sqrt(dx * dx + dy * dy);
    }

    /// @brief Calculates the element stiffness
    /// @param e index of element (rod) in 0 <= e < number_of_elements
    void calculate_kk_element(size_t e);

    /// @brief Calculates the global stiffness
    void calculate_kk();

    /// @brief K matrix for prescribed displacements
    /// @param h is the step-size
    void modify_kk(double h = 1);

    /// @brief CalcDFint calculates internal forces.
    void calculate_delta_ff_int();

    /// @brief Solve solves equilibrium via FEM.
    /// @param number_of_increments number of increments
    void solve(size_t number_of_increments = 1);
};
