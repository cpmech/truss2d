#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../util/doctest.h"

#include "../check/check.h"
#include "../util/printing.h"
#include "truss2d.h"
#include <vector>

using namespace std;

TEST_CASE("truss2d") {

    // input data
    const int nn = 3; // num nodes
    const int nr = 3; // num rods
    double nodes[nn * 2] = {0, 0, 10, 0, 10, 10};
    int con[nr * 2] = {0, 1, 1, 2, 2, 0};
    double props[nr] = {10, 5, 20};
    bool ep[nn * 2] = {true, true, false, true, false, false};
    double ebc[nn * 2] = {0, -0.5, 0, 0.4, 0, 0};
    double nbc[nn * 2] = {0, 0, 0, 0, 2, 1};

    Truss2D truss(nn, nr, nodes, con, ep, nbc, ebc, props);

    truss.Solve();

    // TODO: CHECK(...)
}
