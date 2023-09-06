#include <iostream>

#include "../src/libtruss2d.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) {
    //                           1.0 ^
    //                               |
    //                               |3
    //                               o----> 2.0
    //                             ,'|
    //                           ,'  |
    //                         ,'    |
    //               EA/L=20 ,'      |
    //                 (3) ,'        | EA/L=5
    //                   ,'          | (2)
    //                 ,'            |
    //               ,'              |
    // y           ,'                |
    // |        1,'        (1)       |2
    // |        o--------------------o
    // |____x   ^       EA/L=10      ^
    //         ###                   O
    //         ###                  ---
    // U = {0,-0.5, 0,0.4, -0.5,0.2};
    // F = {-2,-2, 0,1, 2,1};

    /*
        const int nn = 3; // num nodes
        const int nr = 3; // num rods
        double nodes[nn * 2] = {0, 0, 10, 0, 10, 10};
        int con[nr * 2] = {0, 1, 1, 2, 2, 0};
        double props[nr] = {10, 5, 20};
        bool ep[nn * 2] = {true, true, false, true, false, false};
        double ebc[nn * 2] = {0, -0.5, 0, 0.4, 0, 0};
        double nbc[nn * 2] = {0, 0, 0, 0, 2, 1};

        Truss2D truss(nn, nr, nodes, con, ep, nbc, ebc, props);
        truss.solve();

        cout << "==============================================" << endl;
        cout << "                      U                      F" << endl;
        cout << "----------------------------------------------" << endl;
        for (int i = 0; i < truss.Ndof(); i++) {
            print_scientific(truss.U()[i], 23, 15);
            print_scientific(truss.F()[i], 23, 15);
            cout << endl;
        }
        cout << "==============================================" << endl;
        cout << endl;
        */
}
