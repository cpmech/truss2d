#pragma once

#ifdef _MSC_VER
#include <mkl.h>
#define INT MKL_INT
#else
#define INT int
#endif // _MSC_VER

#include <Eigen/Dense>
#include <cmath> // for pow, sqrt, etc
#include <memory>

/// @brief Implements a finite element solver for trusses in 2D
struct Truss2D {
    /// @brief Coordinates x0 y0  x1 y1  ...  xnn ynn (size=2*number_of_nodes)
    Eigen::VectorXd coordinates;

    /// @brief Connectivity 0 1  0 2  1 2  (size=2*number_of_elements)
    Eigen::VectorXi connectivity;

    /// @brief Properties = E*A (size=number_of_elements)
    Eigen::VectorXd properties;

    /// @brief Essential (displacement) prescribed?
    std::vector<bool> essential_prescribed;

    /// @brief Essential (displacement) boundary conditions (size=TOTAL_NDOF)
    Eigen::VectorXd essential_boundary_conditions;

    /// @brief Natural (force) boundary conditions (size=TOTAL_NDOF)
    Eigen::VectorXd natural_boundary_conditions;

    /// @brief Divide E*A by each length?
    bool divide_property_by_length;

    /// @brief Element stiffness matrix (4 x 4)
    Eigen::MatrixXd kk_element;

    /// @brief Global displacements (size=TOTAL_NDOF)
    Eigen::VectorXd uu;

    /// @brief Global forces (size=TOTAL_NDOF)
    Eigen::VectorXd ff;

    /// @brief Internal forces of an element (size=TOTAL_NDOF)
    Eigen::VectorXd ff_int;

    /// @brief Residual = F - F_int (size=TOTAL_NDOF)
    Eigen::VectorXd residual;

    /// @brief Global displacements increments (size=TOTAL_NDOF)
    Eigen::VectorXd delta_uu;

    /// @brief Global forces increments (size=TOTAL_NDOF)
    Eigen::VectorXd delta_ff;

    /// @brief Internal forces increments of an element (size=TOTAL_NDOF)
    Eigen::VectorXd delta_ff_int;

    /// @brief Global stiffness matrix (size=TOTAL_NDOF*TOTAL_NDOF)
    Eigen::MatrixXd kk;

    /// @brief Global stiffness matrix (size=TOTAL_NDOF*TOTAL_NDOF)
    Eigen::MatrixXd kk_copy;

    /// @brief Allocates a new Truss2D structure
    /// @param coordinates x0 y0  x1 y1  ...  xnn ynn (size=2*number_of_nodes)
    /// @param connectivity 0 1  0 2  1 2  (size=2*number_of_elements)
    /// @param properties E*A (size=number_of_elements)
    /// @param essential_prescribed prescribed displacement flags
    /// @param essential_boundary_conditions prescribed displacement values
    /// @param natural_boundary_conditions prescribed external force values
    /// @param divide_property_by_length divide E*A by each rod length
    inline static std::unique_ptr<Truss2D> make_new(
        Eigen::VectorXd coordinates,
        Eigen::VectorXi connectivity,
        Eigen::VectorXd properties,
        std::vector<bool> essential_prescribed,
        Eigen::VectorXd essential_boundary_conditions,
        Eigen::VectorXd natural_boundary_conditions,
        bool divide_property_by_length = true) {
        size_t number_of_nodes = static_cast<size_t>(coordinates.size()) / 2;
        size_t total_ndof = 2 * number_of_nodes;
        return std::unique_ptr<Truss2D>{
            new Truss2D{
                coordinates,
                connectivity,
                properties,
                essential_prescribed,
                essential_boundary_conditions,
                natural_boundary_conditions,
                divide_property_by_length,
                Eigen::MatrixXd(4, 4),                   // kk_element
                Eigen::VectorXd(total_ndof),             // uu
                Eigen::VectorXd(total_ndof),             // ff
                Eigen::VectorXd(total_ndof),             // ff_int
                Eigen::VectorXd(total_ndof),             // residual
                Eigen::VectorXd(total_ndof),             // delta_uu
                Eigen::VectorXd(total_ndof),             // delta_ff
                Eigen::VectorXd(total_ndof),             // delta_ff_int
                Eigen::MatrixXd(total_ndof, total_ndof), // kk
                Eigen::MatrixXd(total_ndof, total_ndof), // kk_copy
            }};
    }

    /// @brief Returns the abscissa of node i (zero based)
    /// @param i index of node in 0 <= i < number_of_nodes
    inline double get_x(size_t i) const {
        if (i >= static_cast<size_t>(coordinates.size())) {
            throw "cannot get x coordinate because the node index is out-of-range";
        }
        return coordinates[i * 2];
    }

    /// @brief Returns the ordinate of node i (zero based)
    /// @param i index of node in 0 <= i < number_of_nodes
    inline double get_y(size_t i) const {
        if (i >= static_cast<size_t>(coordinates.size())) {
            throw "cannot get y coordinate because the node index is out-of-range";
        }
        return coordinates[i * 2 + 1];
    }

    /// @brief Returns the distance between node i and j
    /// @param i index of node in 0 <= i < number_of_nodes
    /// @param j index of node in 0 <= j < number_of_nodes
    inline double calculate_length(size_t i, size_t j) const {
        if (i >= static_cast<size_t>(coordinates.size())) {
            throw "cannot calculate the length because the i-th node index is out-of-range";
        }
        if (j >= static_cast<size_t>(coordinates.size())) {
            throw "cannot calculate the length because the j-th node index is out-of-range";
        }
        return sqrt(pow(get_x(i) - get_x(j), 2.0) + pow(get_y(i) - get_y(j), 2.0));
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
