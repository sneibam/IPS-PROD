#ifndef DEF_DENSITY
#define DEF_DENSITY
/**
 * @file Density.h
 *
 * Interface de la classe Density
 */
#include <iostream>
#include <armadillo>
#include <math.h>
#include <cmath>
#include "Poly.h"
#include "Miscellaneous.h"
#include "Basis.h"

class Density
{
private:
    Basis basis;

public:
    Density();
    //Density(Basis);
    arma::mat calcDensity0(arma::vec, arma::vec);
    arma::mat calcDensityOptimizedv1(arma::vec, arma::vec);
    arma::mat calcDensity1(arma::vec, arma::vec);
    arma::mat calcDensity2(arma::vec, arma::vec);
    arma::mat calcDensity3(arma::vec, arma::vec);
    arma::mat calcDensity4(arma::vec, arma::vec);
    arma::mat calcDensity5(arma::vec, arma::vec);
    arma::cube transformIntoCylindricCoordinates(arma::vec, arma::vec, arma::vec, arma::mat);
};

#endif