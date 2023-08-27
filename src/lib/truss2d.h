#pragma once

#ifdef _MSC_VER
#include <mkl.h>
#define INT MKL_INT
#else
#define INT int
#endif // _MSC_VER

#include <cmath> // for pow, sqrt, etc
#include <memory>

// Implements a Finite Element solver for trusses in 2D
//
// Example of input:
//                             1.0 ^
//                                 |
//                                 |3
//                                 o----> 2.0
//                               ,'|
//                             ,'  |
//                           ,'    |
//                 EA/L=20 ,'      |
//                   (3) ,'        | EA/L=5
//                     ,'          | (2)
//                   ,'            |
//                 ,'              |
//   y           ,'                |
//   |        1,'        (1)       |2
//   |        o--------------------o
//   |____x   ^       EA/L=10      ^
//           ###                   O
//           ###                  ---
//
// const int nn = 3; // num nodes
// const int nr = 3; // num rods
// double nodes     [nn*2] = {0,0, 10,0, 10,10};
// int    con       [nr*2] = {0,1, 1,2, 2,0};
// double props     [nr  ] = {10,5,20};
// bool   ep        [nn*2] = {true,true, false,true, false,false};
// double ebc       [nn*2] = {0,-0.5, 0,0.4, 0,0};
// double nbc       [nn*2] = {0,0, 0,0, 2,1};
// double u_correct [nn*2] = {0,-0.5, 0,0.4, -0.5,0.2};
// double f_correct [nn*2] = {-2,-2, 0,1, 2,1};
//
// Truss2D truss(nn,nr,nodes,con,ep,nbc,ebc,props);
// truss.Solve();
//
// Notation:
//   X,Y : coordinates
//   R   : radii
//   C   : list of connectivity (l:left node, r: right node)
//   P   : properties
//   ep  : is essential prescribed ?
//   bc  : boundary conditions
//   ebc : essential boundary conditions == displacements
//   nbc : natural boundary conditions == forces

/// @brief Implements a finite element solver for trusses in 2D
struct Truss2D {
    int _nn;              // Number of nodes
    int _nr;              // Number of rods
    double const *_nodes; // Coordinates x0 y0  x1 y1  ...  xnn ynn (size=2*nn)
    int const *_con;      // Connectivity 0 1  0 2  1 2  (size=2*nr)
    double const *_props; // Properties (size=nr)
    bool const *_ep;      // Essential (displacement) prescribed?
    double const *_ebc;   // Essential (displacement) boundary conditions
    double const *_nbc;   // Natural (force) boundary conditions
    bool _PbyL;           // Divide P (properties) by L (length) ?
    double _Ke[16];       // Space to hold a stiffness matrix of an element. Set by Ke() method. (must be in col-major format for LAPACK)
    INT _ndof;            // Number of degrees of freedom == 2*nn
    double *_U;           // Global displacements (size=_ndof)
    double *_F;           // Global forces (size=_ndof)
    double *_Fint;        // Internal forces of an element (size=_ndof)
    double *_Res;         // Residual = F - Fint (size=_ndof)
    double *_dU;          // Global displacements increments (size=_ndof)
    double *_dF;          // Global forces increments (size=_ndof)
    double *_dFint;       // Internal forces increments of an element (size=_ndof)
    double *_K;           // Global stiffness matrix. (must be in col-major forces for LAPACK). (size=_ndof*_ndof)
    double *_Kcpy;        // Global stiffness matrix. (must be in col-major forces for LAPACK). (size=_ndof*_ndof)
    INT *_ipiv;           // Index pivot for LAPACK (size=_ndof)

    // Initialize clears U, F, Fint and Res (displacements, forces, internal // forces, and residuals).
    void Initialize();

    // SetNRods (re)set number of rods (elements) (NOTE: this will only work if // Props==NULL).
    void SetNRods(int NRods);

    // SetCon (re)set connectivity (NOTE: this will only work if Props==NULL).
    void SetCon(int const *Con);

    // Destructor
    ~Truss2D();

    // X returns the abscissa of node i
    double X(int i) const { return _nodes[i * 2]; }

    // Y returns the ordinate of node i
    double Y(int i) const { return _nodes[i * 2 + 1]; }

    // L returns the distance between node i and j
    double L(int i, int j) const { return sqrt(pow(X(i) - X(j), 2.0) + pow(Y(i) - Y(j), 2.0)); }

    // Ndof returns the number of degrees of freedom
    int Ndof() const { return _ndof; }

    // CalcKe calculate the element stiffness (Ke).
    //   Input:
    //     e -- index of an element
    void CalcKe(int e);

    // CalcK calculates the global stiffness (K).
    void CalcK();

    // ModifyK modifies K matrix for prescribed displacements.
    //   Input:
    //     h -- step-size
    void ModifyK(double h = 1);

    // CalcDFint calculates internal forces.
    void CalcDFint();

    // Solve solves equilibrium via FEM.
    //   Input:
    //     nInc -- number of increments
    void Solve(int nInc = 1);

    // U returns the displacements
    double const *U() const { return _U; }

    // F returns the external forces
    double const *F() const { return _F; }

    // Fint returns the internal forces
    double const *Fint() const { return _Fint; }

    // Res returns the residual
    double const *Res() const { return _Res; }
};
