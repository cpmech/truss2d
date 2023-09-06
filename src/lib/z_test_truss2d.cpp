#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <vector>

#include "../util/doctest.h"
#include "laclib.h"
#include "truss2d.h"

using namespace std;

TEST_CASE("truss2d") {

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
    // node 0: x-fixed with a vertical
    //         displacement: uy = -0.5
    // node 1: x-fixed with a vertical
    //         displacement: uy = 0.4
    // node 2: fx = 2 and fy = 1
    //
    // EXPECTED RESULTS
    //
    // kk * uu = ff
    //
    // uu_exp = [0.0, -0.5, 0.0, 0.4, -0.5, 0.2]
    // ff_exp = [-2.0, -2.0, 0.0, 1.0, 2.0, 1.0]
}
