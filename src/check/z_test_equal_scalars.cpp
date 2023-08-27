#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../util/doctest.h"

#include "check.h"
#include <string>
#include <vector>

using namespace std;

TEST_CASE("equal_scalars_tol") {
    SUBCASE("float") {
        float a = 1.0;
        float b = 1.0;
        float c = 1.0001;
        float tolerance = 1e-6;
        CHECK(equal_scalars_tol(a, b, tolerance) == true);
        CHECK(equal_scalars_tol(b, c, tolerance) == false);
    }

    SUBCASE("double") {
        double a = 1.0;
        double b = 1.0;
        double c = 1.0000001;
        double tolerance = 1e-15;
        CHECK(equal_scalars_tol(a, b, tolerance) == true);
        CHECK(equal_scalars_tol(b, c, tolerance) == false);
    }
}
