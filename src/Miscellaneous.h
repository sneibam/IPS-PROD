#ifndef DEF_MISCELLANEOUS
#define DEF_MISCELLANEOUS
/**
 * @file Miscellaneous.h
 *
 * Header of the Miscellaneous Class
 *
 */
#include <iostream>
#include <armadillo>
#include <string>
using namespace std;

/**
 * @class Miscellaneous
 *
 * \#Miscellaneous functions for calculus
 */
class Miscellaneous
{
	public:
		static int factorial(int n);
		static arma::mat deriv(arma::mat Y1, arma::mat Y2, arma::mat X1, arma::mat X2);
		static void writeToTxt(arma::mat, std::string);
		static void writeDataSet(arma::vec, arma::vec, arma::mat);
		static std::string cubeToDf3(const arma::cube&);
};

#endif