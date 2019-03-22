
#ifndef DEF_BASIS
#define DEF_BASIS
/**
 * @file Basis.h
 *
 * Interface de la classe Basis
 */
#include <iostream>
#include <armadillo>
#include <math.h>
#include <cmath>
#include "Poly.h"
#include "Miscellaneous.h"


class Basis
{
private:
    double br, bz;
    double calc_nzMax(int N, double Q, int i);
public:
    int mMax;
    arma::ivec nMax;
    arma::imat n_zMax;
    Basis();
    Basis(double, double, int, double);
    arma::vec zPart(arma::vec,int);
    arma::vec rPart(arma::vec, int, int);
    arma::mat basisFunc(int, int, int, arma::vec, arma::vec);
};

#endif