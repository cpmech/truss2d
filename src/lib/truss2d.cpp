#include "truss2d.h"

// LAPACK
extern "C" {
// double-general-solver
#ifdef _MSC_VER
#include <mkl.h>
#else
void dgesv_(int *Np, int *NRHSp, double *A, int *LDAp, int *IPIVp, double *B, int *LDBp, int *INFOp);
#endif
}

void Truss2D::calculate_kk_element(size_t e) {
    size_t number_of_elements = static_cast<size_t>(connectivity.size()) / 2;
    if (e >= number_of_elements) {
        throw "cannot calculate element stiffness because the element index is out-of-range";
    }
    size_t i = connectivity(e * 2);
    size_t j = connectivity(e * 2 + 1);
    double d = calculate_length(i, j);
    double c = (get_x(j) - get_x(i)) / d;
    double s = (get_y(j) - get_y(i)) / d;
    double p = (divide_property_by_length ? properties[e] / d : properties[e]);

    kk_element(0, 0) = p * c * c;
    kk_element(0, 1) = p * c * s;
    kk_element(0, 2) = -p * c * c;
    kk_element(0, 3) = -p * c * s;

    kk_element(1, 0) = kk_element(0, 1);
    kk_element(1, 1) = p * s * s;
    kk_element(1, 2) = kk_element(0, 3);
    kk_element(1, 3) = -p * s * s;

    kk_element(2, 0) = kk_element(0, 2);
    kk_element(2, 1) = kk_element(0, 3);
    kk_element(2, 2) = kk_element(0, 0);
    kk_element(2, 3) = kk_element(0, 1);

    kk_element(3, 0) = kk_element(0, 3);
    kk_element(3, 1) = kk_element(1, 3);
    kk_element(3, 2) = kk_element(0, 1);
    kk_element(3, 3) = kk_element(1, 1);
}

void Truss2D::calculate_kk() {
    size_t number_of_nodes = static_cast<size_t>(coordinates.size()) / 2;
    size_t total_ndof = 2 * number_of_nodes;
    for (size_t i = 0; i < total_ndof; ++i) {
        for (size_t j = 0; j < total_ndof; ++j) {
            kk(i, j) = 0.0;
        }
    }
    size_t number_of_elements = static_cast<size_t>(connectivity.size()) / 2;
    for (size_t e = 0; e < number_of_elements; ++e) {
        calculate_kk_element(e);
        size_t l = connectivity(e * 2);                     // left node
        size_t r = connectivity(e * 2 + 1);                 // right node
        size_t m[4] = {l * 2, l * 2 + 1, r * 2, r * 2 + 1}; // map (local=>global)
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                kk(m[i], m[j]) += kk_element(i, j);
            }
        }
    }
}

void Truss2D::modify_kk(double h) {
    size_t number_of_nodes = static_cast<size_t>(coordinates.size()) / 2;
    size_t total_ndof = 2 * number_of_nodes;

    // set prescribed values
    for (size_t i = 0; i < total_ndof; ++i) {
        delta_ff(i) = 0.0;
        delta_uu(i) = 0.0;
        delta_ff_int[i] = 0.0;
        if (essential_prescribed[i]) {
            delta_uu(i) = essential_boundary_conditions(i) * h;
        } else {
            delta_ff(i) = natural_boundary_conditions(i) * h;
            delta_uu(i) = delta_ff[i];
            for (size_t j = 0; j < total_ndof; ++j) {
                if (essential_prescribed[j]) {
                    delta_uu(i) -= kk(i, j) * essential_boundary_conditions(j) * h;
                }
            }
        }
    }

    // clear lines and columns of K for prescribed displacements => modified stiffness
    for (size_t i = 0; i < total_ndof; ++i) {
        if (essential_prescribed[i]) {
            for (size_t j = 0; j < total_ndof; ++j) {
                kk(i, j) = 0.0;
                kk(j, i) = 0.0;
            }
            kk(i, i) = 1.0;
        }
    }
}

void Truss2D::calculate_delta_ff_int() {
    size_t number_of_elements = static_cast<size_t>(connectivity.size()) / 2;
    for (size_t e = 0; e < number_of_elements; ++e) {
        size_t i = connectivity(e * 2);
        size_t j = connectivity(e * 2 + 1);
        double d = calculate_length(i, j);
        double c = (get_x(j) - get_x(i)) / d;
        double s = (get_y(j) - get_y(i)) / d;
        double p = (divide_property_by_length ? properties[e] / d : properties[e]);
        double delta_length = c * delta_uu[j * 2] + s * delta_uu[j * 2 + 1] - (c * delta_uu[i * 2] + s * delta_uu[i * 2 + 1]);
        double delta_axial_force = delta_length * p;
        delta_ff_int(i * 2) += -c * delta_axial_force;
        delta_ff_int(i * 2 + 1) += -s * delta_axial_force;
        delta_ff_int(j * 2) += c * delta_axial_force;
        delta_ff_int(j * 2 + 1) += s * delta_axial_force;
    }
}

void Truss2D::solve(size_t number_of_increments) {
    size_t number_of_nodes = static_cast<size_t>(coordinates.size()) / 2;
    size_t total_ndof = 2 * number_of_nodes;
    double h = 1.0 / number_of_increments;
    for (size_t i = 0; i < number_of_increments; ++i) {
        // assembly
        calculate_kk();

        // save a copy of K for later recovering dF (could save only K21 and K22)
        for (size_t j = 0; j < total_ndof; ++j) {
            for (size_t i = 0; i < total_ndof; ++i) {
                kk_copy(i, j) = kk(i, j);
            }
        }

        // modify K for prescribed displacements
        modify_kk(h);

        // solve dU = inv(K)*dF
        // INT info = 0;
        // INT nrhs = 1;
        // dgesv_(&_ndof, &nrhs, _K, &_ndof, _ipiv, _dU, &_ndof, &info);

        // solve external forces increments for prescribed essential (displacement) dofs
        for (size_t j = 0; j < total_ndof; ++j) {
            if (essential_prescribed[j]) {
                for (size_t m = 0; m < total_ndof; ++m) {
                    delta_ff[j] += kk_copy(j, m) * delta_uu(m);
                }
            }
        }

        // internal forces
        calculate_delta_ff_int();

        // update
        for (size_t j = 0; j < total_ndof; ++j) {
            uu(j) += delta_uu(j);
            ff(j) += delta_ff(j);
            ff_int(j) += delta_ff_int(j);
            residual(j) = ff(j) - ff_int(j);
        }
    }
}
