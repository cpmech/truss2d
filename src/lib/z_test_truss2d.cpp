#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../util/doctest.h"

#include "../check/check.h"
#include "../util/printing.h"
#include "truss2d.h"
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

TEST_CASE("truss2d") {

    size_t number_of_nodes = 3;
    size_t number_of_elements = 3;

    VectorXd coordinates(2 * number_of_nodes);
    coordinates << 0.0, 0.0, 10.0, 0.0, 10.0, 10.0;

    VectorXi connectivity(2 * number_of_elements);
    connectivity << 0, 1, 1, 2, 2, 0;

    VectorXd properties(number_of_elements);
    properties << 10.0, 5.0, 20.0;

    vector<bool> essential_prescribed = {true, true, false, true, false, false};

    VectorXd essential_boundary_conditions(2 * number_of_nodes);
    essential_boundary_conditions << 0.0, -0.5, 0.0, 0.4, 0.0, 0.0;

    VectorXd natural_boundary_conditions(2 * number_of_nodes);
    natural_boundary_conditions << 0.0, 0.0, 0.0, 0.0, 2.0, 1.0;

    auto truss = Truss2D::make_new(
        coordinates,
        connectivity,
        properties,
        essential_prescribed,
        essential_boundary_conditions,
        natural_boundary_conditions);

    truss->solve();

    // TODO: CHECK(...)
}
