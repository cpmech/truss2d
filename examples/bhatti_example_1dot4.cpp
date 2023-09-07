#include <iostream>

#include "../src/libtruss2d.h"
#include "laclib.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) {

    // Bhatti's Example 1.4 on page 25
    //
    // Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
    //
    // TEST GOAL
    //
    // This test verifies a 2D frame with rod elements and concentrated forces
    //
    // MESH
    //
    //               (3)
    //               [2]
    //     2----------------------3
    //     |'.  (4)           _.-'
    //     |  '.[3]       _.-'
    //     |    '.    _.-'  (1)
    // (2) |      '1-'      [1]
    // [2] |      /
    //     |     /
    //     |    / (0)   The lines are ROD (Lin2) elements
    //     |   /  [1]
    //     |  /
    //     | /    (#) indicates cell id
    //     0'     [#] indicates attribute id
    //
    // BOUNDARY CONDITIONS
    //
    // Fully fixed @ points 0 and 3
    // Concentrated load @ point 1 with Fy = -150,000
    //
    // CONFIGURATION AND PARAMETERS
    //
    // Static simulation
    // Attribute 1: Area = 4,000; Young = 200,000
    // Attribute 2: Area = 3,000; Young = 200,000
    // Attribute 3: Area = 2,000; Young =  70,000

    // nodes
    auto coordinates = vector<double>{
        0.0, 0.0,    // 0
        1500, 3500,  // 1
        0.0, 5000,   // 2
        5000, 5000}; // 3

    // bars
    auto connectivity = vector<size_t>{
        0, 1,  // 0
        1, 3,  // 1
        0, 2,  // 2
        2, 3,  // 3
        1, 2}; // 4

    // E*A
    double mat1 = 200000.0 * 4000.0;
    double mat2 = 200000.0 * 3000.0;
    double mat3 = 70000.0 * 2000.0;
    auto properties = vector<double>{mat1, mat1, mat2, mat2, mat3};

    // boundary conditions
    map<node_dof_pair_t, double> essential_bcs{
        {{0, AlongX}, 0.0},
        {{0, AlongY}, 0.0},
        {{3, AlongX}, 0.0},
        {{3, AlongY}, 0.0}};
    map<node_dof_pair_t, double> natural_bcs{
        {{1, AlongY}, -150000}};

    // bhatti's solution
    auto correct_uu = vector<double>{
        0.000000000000000e00, 0.000000000000000e00,    // 0
        5.389536380057676e-01, -9.530613006371175e-01, // 1
        2.647036149579491e-01, -2.647036149579490e-01, // 2
        0.000000000000000e00, 0.000000000000000e00};   // 3

    // solve
    auto truss = Truss2d::make_new(coordinates, connectivity, properties, essential_bcs, natural_bcs);
    truss->solve();
    print_vector("uu", truss->uu);

    // check
    if (equal_vectors_tol(truss->uu, correct_uu, 1e-15)) {
        cout << "OK" << endl;
        return 0;
    } else {
        cout << "FAIL" << endl;
        return 1;
    }
}
