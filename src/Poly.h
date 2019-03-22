/**
 * @file Calcul.h
 *
 * Header of Calcul class
 *
 */

#ifndef DEF_POLY
#define DEF_POLY
#include <armadillo>
#include <vector>
using namespace arma;
using namespace std;

/**
 * @class Poly
 *
 * This Class computes the Hermite Polynomial and the Wave function
 *
 */

class Poly
{
	public:
		/// Default constructor of the object
		Poly();

		/// This function computes the Hermite Polynomial using this recurrence relation
		void calcHermite(int,vec);
		void calcLaguerre(int,int,vec);
		vec hermite(int);
		vec laguerre(int,int);

		/// Destructor of the object
		~Poly();

	private:
		mat hermiteVals	;
		cube laguerreVals;
};

#endif
