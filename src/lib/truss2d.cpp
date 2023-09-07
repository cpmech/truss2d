#include "truss2d.h"

void Truss2d::calculate_element_stiffness(size_t e) {
    if (e >= number_of_elements) {
        throw "cannot calculate element stiffness because the element index is out-of-range";
    }
    size_t a = connectivity[e * 2];
    size_t b = connectivity[e * 2 + 1];
    double l = calculate_length(a, b);
    auto xa = coordinates[a * 2];
    auto ya = coordinates[a * 2 + 1];
    auto xb = coordinates[b * 2];
    auto yb = coordinates[b * 2 + 1];
    double c = (xb - xa) / l;
    double s = (yb - ya) / l;
    double p = properties[e] / l;

    // computing upper triangle only
    //      _                   _
    //     |  c*c c*s -c*c -c*s  | 0
    // E A |   .  s*s -c*s -s*s  | 1
    // --- |   .   .   c*c  c*s  | 2
    //  L  |_  .   .    .   s*s _| 3
    //         0   1    2    3
    kk_element->set(0, 0, p * c * c);
    kk_element->set(0, 1, p * c * s);
    kk_element->set(0, 2, -p * c * c);
    kk_element->set(0, 3, -p * c * s);

    kk_element->set(1, 1, p * s * s);
    kk_element->set(1, 2, -p * c * s);
    kk_element->set(1, 3, -p * s * s);

    kk_element->set(2, 2, p * c * c);
    kk_element->set(2, 3, p * c * s);

    kk_element->set(3, 3, p * s * s);
}

void Truss2d::calculate_rhs_and_global_stiffness() {
    // initialize uu and right-hand side vector
    for (size_t i = 0; i < total_ndof; ++i) {
        if (essential_prescribed[i]) {
            uu[i] = essential_boundary_conditions[i];  // needed to correct RHS vector
            rhs[i] = essential_boundary_conditions[i]; // because diagonal(K;prescribed) = 1
        } else {
            uu[i] = 0.0;                             // irrelevant, actually
            rhs[i] = natural_boundary_conditions[i]; // external forces
        }
    }

    // assemble stiffness and fix RHS vector
    kk_coo->pos = 0; // reset
    for (size_t e = 0; e < number_of_elements; ++e) {
        calculate_element_stiffness(e);
        size_t a = connectivity[e * 2];
        size_t b = connectivity[e * 2 + 1];
        size_t m[4] = {a * 2, a * 2 + 1, b * 2, b * 2 + 1}; // map local => global
        for (size_t i = 0; i < 4; ++i) {
            if (!essential_prescribed[m[i]]) {
                // correct RHS
                for (size_t j = 0; j < 4; ++j) {
                    if (essential_prescribed[m[j]]) {
                        if (j >= i) {
                            rhs[m[i]] -= kk_element->get(i, j) * uu[m[j]];
                        } else {
                            // must get (i,j) from upper triangle
                            rhs[m[i]] -= kk_element->get(j, i) * uu[m[j]];
                        }
                    }
                }
                // assemble upper triangle into global stiffness
                for (size_t j = i; j < 4; ++j) { // j = i => local upper triangle
                    if (!essential_prescribed[m[j]]) {
                        if (m[j] >= m[i]) {
                            kk_coo->put(m[i], m[j], kk_element->get(i, j));
                        } else {
                            // must go to the global upper triangle
                            kk_coo->put(m[j], m[i], kk_element->get(i, j));
                        }
                    }
                }
            }
        }
    }

    // put ones on the diagonal of the global stiffness matrix
    for (size_t i = 0; i < total_ndof; i++) {
        if (essential_prescribed[i]) {
            kk_coo->put(i, i, 1.0);
        }
    }

    // convert COO to CSR
    if (kk_csr == NULL) {
        kk_csr = CsrMatrixMkl::from(kk_coo);
    } else {
        kk_csr.reset();
        kk_csr = CsrMatrixMkl::from(kk_coo);
    }
}

void Truss2d::solve() {
    lin_sys_solver->analyze(kk_csr);
    lin_sys_solver->factorize(kk_csr);
    lin_sys_solver->solve(uu, rhs); // uu = inv(kk) * ff
}
