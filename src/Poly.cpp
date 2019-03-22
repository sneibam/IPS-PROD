/**
 * @file Calcul.cpp
 *
 * Implementation of the class Calcul
 *
 */
#include <stdlib.h>
#include "Poly.h"
#include <math.h>
#include <cmath>


/** Default Constructor of Calcul class
 *
 *
 */
Poly::Poly()
{

}

/** Computing the Hermite Polynomial
 *
 *
 * This function computes the Hermite Polynomial using the reccurence formula
 *
 * \#Definition(physicists version):
 * \f[\forall n\ge 0, H_n(z)\equiv (-1)^n e^{z^2} \frac{d^n}{dz^n}\left( e^{-z^2} \right).\f]
 *
 * \#Reccurence relation:
 * \f[H_0(z) = 1\f]
 * \f[H_1(z) = 2z\f]
 * \f[\forall n\ge 1, H_{n+1}(z) = 2zH_n(z)-2nH_{n-1}(z).\f]
 *
 * @return a matrix containing the values of the function
 */


void Poly::calcHermite(int n_max, vec Z)
{

    int i;

    if (n_max < 0)
    {
        std::cout << "Error !" << std::endl;
        exit(0);
    }

    hermiteVals.ones(Z.n_rows,n_max+1);


    if (n_max == 0)
    {
        return;
    }

    hermiteVals.col(1) = 2 * Z;

    for (i = 2; i <= n_max; i++)
    {
        hermiteVals.col(i) = (2 * Z % hermiteVals.col(i - 1)) - (2 * (i-1) * hermiteVals.col(i - 2));
    }


}

void Poly::calcLaguerre(int m, int n, vec Z) 
{

	int i;
    laguerreVals.zeros(Z.n_rows,m+1,n+1);

    //cas n=0
    laguerreVals.slice(0) = laguerreVals.slice(0) + 1;

    //cas n=1
    mat L1 = zeros(Z.n_rows, m+1);
    mat matZ = zeros(size(L1));
    mat M = zeros(size(L1));
    for (i=0; i<=m; i++)
    {
        L1.col(i) = 1 + i - Z;
        matZ.col(i) = Z;
        M.col(i).fill(i);
    }
    laguerreVals.slice(1) = L1;

    //les autres cas
    for (i=2; i<=n; i++)
    {
		laguerreVals.slice(i) = (2 + (M - 1 - matZ) / i) % laguerreVals.slice(i-1) - (1 + (M - 1)/ i) % laguerreVals.slice(i-2);
	}
}

vec Poly::hermite(int n)
{
	return hermiteVals.col(n);
}

vec Poly::laguerre(int m, int n)
{
	return laguerreVals.slice(n).col(m);
}

/**
*
* Destructor of the object Calcul
*/

Poly::~Poly()
{

}
