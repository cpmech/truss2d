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

/*
Truss2D::Truss2D(int NNodes, int NRods, double const *Nodes, int const *Connects, bool const *EssenPresc,
                 double const *NaturalBC, double const *EssenBC, double const *Props, bool DividePbyL) {
    // set input data
    _nn = NNodes;
    _nr = NRods;
    _nodes = Nodes;
    _con = Connects;
    _props = Props;
    _ep = EssenPresc;
    _ebc = EssenBC;
    _nbc = NaturalBC;
    _PbyL = DividePbyL;

    // allocate arrays
    _ndof = 2 * _nn;
    _U = new double[_ndof];
    _F = new double[_ndof];
    _Fint = new double[_ndof];
    _Res = new double[_ndof];
    _dU = new double[_ndof];
    _dF = new double[_ndof];
    _dFint = new double[_ndof];
    _K = new double[_ndof * _ndof];
    _Kcpy = new double[_ndof * _ndof];
    _ipiv = new INT[_ndof];

    // initialize arrays
    Initialize();
}
*/

void Truss2D::Initialize() {
    for (int i = 0; i < _ndof; ++i) {
        _U[i] = 0.0;
        _F[i] = 0.0;
        _Fint[i] = 0.0;
        _Res[i] = 0.0;
    }
}

void Truss2D::SetNRods(int NRods) { _nr = NRods; }

void Truss2D::SetCon(int const *Con) { _con = Con; }

Truss2D::~Truss2D() {
    delete[] _U;
    delete[] _F;
    delete[] _Fint;
    delete[] _Res;
    delete[] _dU;
    delete[] _dF;
    delete[] _dFint;
    delete[] _K;
    delete[] _Kcpy;
    delete[] _ipiv;
}

void Truss2D::CalcKe(int e) {
    int i = _con[e * 2];
    int j = _con[e * 2 + 1];
    double d = L(i, j);
    double c = (X(j) - X(i)) / d;
    double s = (Y(j) - Y(i)) / d;
    double p = (_props == NULL ? 1.0 : (_PbyL ? _props[e] / d : _props[e]));
    _Ke[0] = p * c * c;
    _Ke[4] = p * c * s;
    _Ke[8] = -p * c * c;
    _Ke[12] = -p * c * s;
    _Ke[1] = _Ke[4];
    _Ke[5] = p * s * s;
    _Ke[9] = _Ke[12];
    _Ke[13] = -p * s * s;
    _Ke[2] = _Ke[8];
    _Ke[6] = _Ke[9];
    _Ke[10] = _Ke[0];
    _Ke[14] = _Ke[4];
    _Ke[3] = _Ke[12];
    _Ke[7] = _Ke[13];
    _Ke[11] = _Ke[14];
    _Ke[15] = _Ke[5];
}

void Truss2D::CalcK() {
    for (int i = 0; i < _ndof * _ndof; ++i) {
        _K[i] = 0.0;
    }
    for (int e = 0; e < _nr; ++e) {
        CalcKe(e);
        int l = _con[e * 2];                             // left node
        int r = _con[e * 2 + 1];                         // right node
        int m[4] = {l * 2, l * 2 + 1, r * 2, r * 2 + 1}; // map (local=>global)
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                _K[m[j] * _ndof + m[i]] += _Ke[j * 4 + i];
    }
}

void Truss2D::ModifyK(double h) {
    // Set dF and dU (workspace) vectors
    for (int i = 0; i < _ndof; ++i) {
        _dF[i] = 0.0;
        _dU[i] = 0.0;
        _dFint[i] = 0.0;
    }
    for (int i = 0; i < _ndof; ++i) {
        if (_ep[i]) {
            _dU[i] = (_ebc == NULL ? 0.0 : _ebc[i] * h);
        } else {
            _dF[i] = _dU[i] = _nbc[i] * h;
            for (int j = 0; j < _ndof; ++j) {
                if (_ep[j]) {
                    _dU[i] -= (_ebc == NULL ? 0.0 : _K[j * _ndof + i] * _ebc[j] * h);
                }
            }
        }
    }

    // Clear lines and columns of K for prescribed displacements => modified
    // stiffness
    for (int i = 0; i < _ndof; ++i) {
        if (_ep[i]) {
            for (int j = 0; j < _ndof; ++j) {
                _K[j * _ndof + i] = 0.0;
                _K[i * _ndof + j] = 0.0;
            }
            _K[i * _ndof + i] = 1.0;
        }
    }
}

void Truss2D::CalcDFint() {
    for (int e = 0; e < _nr; ++e) {
        int i = _con[e * 2];
        int j = _con[e * 2 + 1];
        double d = L(i, j);
        double c = (X(j) - X(i)) / d;
        double s = (Y(j) - Y(i)) / d;
        double p = (_props == NULL ? 1.0 : (_PbyL ? _props[e] / d : _props[e]));
        double elong = c * _dU[j * 2] + s * _dU[j * 2 + 1] - (c * _dU[i * 2] + s * _dU[i * 2 + 1]); // elongation
        double dfn = elong * p;                                                                     // normal force increment
        _dFint[i * 2] += -c * dfn;
        _dFint[i * 2 + 1] += -s * dfn;
        _dFint[j * 2] += c * dfn;
        _dFint[j * 2 + 1] += s * dfn; // global coordinates
    }
}

void Truss2D::Solve(int nInc) {
    double h = 1.0 / nInc;
    for (int i = 0; i < nInc; ++i) {
        // Assembly
        CalcK();

        // Save a copy of K for later recovering of dF (could save only K21 and K22)
        for (int j = 0; j < _ndof * _ndof; ++j) {
            _Kcpy[j] = _K[j];
        }

        // Modify K for prescribed displacements
        ModifyK(h);

        // Solve dU = inv(K)*dF
        INT info = 0;
        INT nrhs = 1;
        dgesv_(&_ndof, &nrhs, _K, &_ndof, _ipiv, _dU, &_ndof, &info);

        // Solve external forces increments for prescribed essential (displacement) dofs
        for (int j = 0; j < _ndof; ++j) {
            if (_ep[j]) {
                for (int m = 0; m < _ndof; ++m) {
                    _dF[j] += _Kcpy[m * _ndof + j] * _dU[m];
                }
            }
        }

        // Calc internal forces
        CalcDFint();

        // Update
        for (int j = 0; j < _ndof; ++j) {
            _U[j] += _dU[j];
            _F[j] += _dF[j];
            _Fint[j] += _dFint[j];
            _Res[j] = _F[j] - _Fint[j];
        }
    }
}
