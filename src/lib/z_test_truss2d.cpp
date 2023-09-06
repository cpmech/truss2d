#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <map>
#include <vector>

#include "../util/doctest.h"
#include "laclib.h"
#include "truss2d.h"

using namespace std;

// sqrt(2) <https://oeis.org/A002193>
#define SQRT_2 1.41421356237309504880168872420969807856967187537694807317667973799073247846210703885038753432764157

TEST_CASE("truss2d") {
    SUBCASE("three-member truss") {
        // GEOMETRY
        //
        //                      fy=1 ↑
        // ---                       2 →
        //  ↑                      ,'| fx=2
        //  |                    ,'  |
        //  |                  ,'    |
        //  |       EA=200√2 ,'      |
        // 10          (2) ,'        | EA=50
        //  |            ,'          | (1)
        //  |          ,'            |
        //  |        ,'              |
        //  |      ,'    EA=100      |
        //  ↓    ,'       (0)        |
        // ---  0--------------------1
        //     | |                  | |
        //      ⇊ uy=-0.5     uy=0.4 ⇈
        //
        //      |←------- 10 -------→|
        //
        // BOUNDARY CONDITIONS
        //
        // node 0: x-fixed with a vertical displacement: uy = -0.5
        // node 1: x-fixed with a vertical displacement: uy = 0.4
        // node 2: fx = 2.0 and fy = 1.0
        //
        // EXPECTED RESULTS
        //
        // kk * uu = ff
        //
        // correct_uu = {0.0, -0.5, 0.0, 0.4, -0.5, 0.2}
        // correct_ff = {-2.0, -2.0, 0.0, 1.0, 2.0, 1.0}
        //
        // Reference
        // Carlos Felippa I-FEM Page 3-12 Chapter 3 The Direct Stiffness Method II

        // input data
        auto coordinates = vector<double>{0.0, 0.0, 10.0, 0.0, 10.0, 10.0};
        auto connectivity = vector<size_t>{0, 1, 1, 2, 2, 0};
        auto properties = vector<double>{100.0, 50.0, 200.0 * SQRT_2};
        map<node_dof_pair_t, double> essential_bcs{
            {{0, AlongX}, 0.0},
            {{0, AlongY}, -0.5},
            {{1, AlongY}, 0.4}};
        map<node_dof_pair_t, double> natural_bcs{
            {{2, AlongX}, 2.0},
            {{2, AlongY}, 1.0},
        };

        // allocate truss solver
        auto truss = Truss2d::make_new(coordinates, connectivity, properties, essential_bcs, natural_bcs);

        // check boundary condition arrays
        auto correct_ep = vector<bool>{true, true, false, true, false, false};
        auto correct_ebc = vector<double>{0.0, -0.5, 0.0, 0.4, 0.0, 0.0};
        auto correct_nbc = vector<double>{0.0, 0.0, 0.0, 0.0, 2.0, 1.0};
        CHECK(equal_vectors(truss->essential_prescribed, correct_ep));
        CHECK(equal_vectors_tol(truss->essential_boundary_conditions, correct_ebc, 1e-17));
        CHECK(equal_vectors_tol(truss->natural_boundary_conditions, correct_nbc, 1e-17));

        // check element length
        CHECK(equal_scalars_tol(truss->calculate_length(0, 1), 10.0, 1e-15));
        CHECK(equal_scalars_tol(truss->calculate_length(0, 2), 10.0 * SQRT_2, 1e-15));
        CHECK(equal_scalars_tol(truss->calculate_length(1, 2), 10.0, 1e-15));

        // check element stiffness
        truss->calculate_element_stiffness(0);
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 0), 10.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 1), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 2), -10.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 1), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 2), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 2), 10.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(3, 3), 0.0, 1e-15));

        truss->calculate_element_stiffness(1);
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 0), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 1), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 2), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 1), 5.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 2), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 3), -5.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 2), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(3, 3), 5.0, 1e-15));

        truss->calculate_element_stiffness(2);
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 0), 10.0, 1e-14));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 1), 10.0, 1e-14));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 2), -10.0, 1e-14));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 3), -10.0, 1e-14));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 1), 10.0, 1e-14));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 2), -10.0, 1e-14));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 3), -10.0, 1e-14));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 2), 10.0, 1e-14));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 3), 10.0, 1e-14));
        CHECK(equal_scalars_tol(truss->kk_element->get(3, 3), 10.0, 1e-14));

        // check global stiffness matrix
        truss->calculate_global_stiffness();
        auto kk = truss->kk_coo->as_matrix();
        // kk->print();
        CHECK(equal_scalars_tol(kk->get(0, 0), 20.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(0, 1), 10.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(0, 2), -10.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(0, 3), 0.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(0, 4), -10.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(0, 5), -10.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(1, 1), 10.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(1, 2), 0.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(1, 3), 0.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(1, 4), -10.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(1, 5), -10.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(2, 2), 10.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(2, 3), 0.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(2, 4), 0.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(2, 5), 0.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(3, 3), 5.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(3, 4), 0.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(3, 5), -5.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(4, 4), 10.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(4, 5), 10.0, 1e-14));
        CHECK(equal_scalars_tol(kk->get(5, 5), 15.0, 1e-14));

        // solve mechanical problem

        // check solution
        auto correct_u = vector<double>{0.0, -0.5, 0.0, 0.4, -0.5, 0.2};
        auto correct_f = vector<double>{-2.0, -2.0, 0.0, 1.0, 2.0, 1.0};
    }

    SUBCASE("eight-member truss") {
        // GEOMETRY
        //
        //                                     (3)
        // ---                       2--------------------4 → fx=50
        //  ↑    E = 30000         ,'|'.                .'|
        //  |    A = 10          ,'  |  '.            .'  |
        //  |                  ,'    |    '.     (5).'    |
        //  |                ,'      |      '.    .'      |
        // 144         (0) ,'        |        '. '        |(7)
        //  |            ,'          |(2)     . '.        |
        //  |          ,'            |      .'    '.      |
        //  |        ,'              |    .'     (4)'.    |
        //  |      ,'                |  .'            '.  |
        //  ↓    ,'       (1)        |.'      (6)       '.|
        // ---  0--------------------1--------------------3
        //    fixed                  ↓ fy=-100          fixed
        //
        //      |←------ 192 -------→|←------ 192 -------→|
        //
        // BOUNDARY CONDITIONS
        //
        // node 0: x and y fixed
        // node 3: x and y fixed
        // node 1: vertical force fy = -100
        // node 4: horizontal force fx = 50
        //
        // EXPECTED RESULTS
        //
        // kk * uu = ff
        //
        // correct_uu = TODO
        // correct_ff = TODO
        //
        // REFERENCE
        // CEE 421L. Matrix Structural Analysis – Duke University – Fall 2014 – H.P. Gavin

        // input data
        auto coordinates = vector<double>{0.0, 0.0, 192.0, 0.0, 192.0, 144.0, 384.0, 0.0, 384.0, 144.0};
        auto connectivity = vector<size_t>{0, 2, 0, 1, 1, 2, 2, 4, 2, 3, 1, 4, 1, 3, 3, 4};
        auto properties = vector<double>(8, 30000.0 * 10.0);
        map<node_dof_pair_t, double> essential_bcs{
            {{0, AlongX}, 0.0},
            {{0, AlongY}, 0.0},
            {{3, AlongX}, 0.0},
            {{3, AlongY}, 0.0}};
        map<node_dof_pair_t, double> natural_bcs{
            {{1, AlongY}, -100.0},
            {{4, AlongX}, 50.0},
        };

        // allocate truss solver
        auto truss = Truss2d::make_new(coordinates, connectivity, properties, essential_bcs, natural_bcs);

        // check element length
        CHECK(equal_scalars_tol(truss->calculate_length(0, 1), 192.0, 1e-15));
        CHECK(equal_scalars_tol(truss->calculate_length(0, 2), 240.0, 1e-15));
        CHECK(equal_scalars_tol(truss->calculate_length(1, 4), 240.0, 1e-15));
        CHECK(equal_scalars_tol(truss->calculate_length(2, 3), 240.0, 1e-15));
        CHECK(equal_scalars_tol(truss->calculate_length(3, 2), 240.0, 1e-15));
        CHECK(equal_scalars_tol(truss->calculate_length(4, 3), 144.0, 1e-15));

        // check element stiffness
        truss->calculate_element_stiffness(0);
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 0), 800.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 1), 600.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 2), -800.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 3), -600.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 0), 0.0, 1e-15)); // not 600, because of using upper triangle only
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 1), 450.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 2), -600.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 3), -450.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 2), 800.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 3), 600.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(3, 3), 450.0, 1e-15));

        truss->calculate_element_stiffness(1);
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 0), 1562.5, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 1), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 2), -1562.5, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 1), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 2), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 2), 1562.5, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(3, 3), 0.0, 1e-15));

        truss->calculate_element_stiffness(2);
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 0), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 1), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 2), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 1), 2083.3333333333333, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 2), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 3), -2083.3333333333333, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 2), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(3, 3), 2083.3333333333333, 1e-15));

        truss->calculate_element_stiffness(3); // same as rod # 1
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 0), 1562.5, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 1), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 2), -1562.5, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 1), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 2), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 2), 1562.5, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(3, 3), 0.0, 1e-15));

        truss->calculate_element_stiffness(4);
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 0), 800.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 1), -600.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 2), -800.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 3), 600.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 1), 450.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 2), 600.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 3), -450.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 2), 800.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 3), -600.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(3, 3), 450.0, 1e-15));

        truss->calculate_element_stiffness(5); // same as rod # 0
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 0), 800.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 1), 600.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 2), -800.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 3), -600.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 1), 450.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 2), -600.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 3), -450.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 2), 800.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 3), 600.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(3, 3), 450.0, 1e-15));

        truss->calculate_element_stiffness(6); // same as rod # 1
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 0), 1562.5, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 1), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 2), -1562.5, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 1), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 2), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 2), 1562.5, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(3, 3), 0.0, 1e-15));

        truss->calculate_element_stiffness(7); // same as rod # 2
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 0), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 1), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 2), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(0, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 1), 2083.3333333333333, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 2), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(1, 3), -2083.3333333333333, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 2), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(2, 3), 0.0, 1e-15));
        CHECK(equal_scalars_tol(truss->kk_element->get(3, 3), 2083.3333333333333, 1e-15));

        // check global stiffness matrix
        truss->calculate_global_stiffness();
        auto kk = truss->kk_coo->as_matrix();
        // kk->print();
        // truss->kk_csr->print("global stiffness matrix");
        CHECK(equal_scalars_tol(kk->get(2, 2), 3925.0, 1e-15));
        CHECK(equal_scalars_tol(kk->get(2, 3), 600.0, 1e-15));
        CHECK(equal_scalars_tol(kk->get(2, 4), 0.0, 1e-15));
        CHECK(equal_scalars_tol(kk->get(2, 5), 0.0, 1e-15));
        CHECK(equal_scalars_tol(kk->get(2, 8), -800.0, 1e-15));
        CHECK(equal_scalars_tol(kk->get(2, 9), -600.0, 1e-15));
        CHECK(equal_scalars_tol(kk->get(3, 3), 2533.3333333333333, 1e-15));
        CHECK(equal_scalars_tol(kk->get(3, 4), 0.0, 1e-15));
        CHECK(equal_scalars_tol(kk->get(3, 5), -2083.3333333333333, 1e-15));
        CHECK(equal_scalars_tol(kk->get(3, 8), -600.0, 1e-15));
        CHECK(equal_scalars_tol(kk->get(3, 9), -450.0, 1e-15));
        CHECK(equal_scalars_tol(kk->get(4, 4), 3162.50, 1e-15));
        CHECK(equal_scalars_tol(kk->get(4, 5), 0.0, 1e-15));
        CHECK(equal_scalars_tol(kk->get(4, 8), -1562.50, 1e-15));
        CHECK(equal_scalars_tol(kk->get(4, 9), 0.0, 1e-15));
        CHECK(equal_scalars_tol(kk->get(5, 5), 2983.3333333333333, 1e-15));
        CHECK(equal_scalars_tol(kk->get(5, 8), 0.0, 1e-15));
        CHECK(equal_scalars_tol(kk->get(5, 9), 0.0, 1e-15));
        CHECK(equal_scalars_tol(kk->get(8, 8), 2362.50, 1e-15));
        CHECK(equal_scalars_tol(kk->get(8, 9), 600.0, 1e-15));
        CHECK(equal_scalars_tol(kk->get(9, 9), 2533.3333333333333, 1e-15));
    }
}
