/**
 * @file Miscellaneous.cpp
 *
 * \#Implementation of Miscellaneous functions for calculus
 */

#include "Miscellaneous.h"


/** Computing the factorial of n
*
* This function computes the factorial of n using reccurence
*
* @param n integer
*
* @return the result of the factorial of n
*/
int Miscellaneous::factorial(int n)
{
	if(n == 0)
		return 1;
	else if(n == 1)
		return 1;
	else if(n == 2)
		return 2;
	else
		return n * factorial(n - 1); // recursive call to factorial()
}


/**
 * Computing the derivative of a function
 *
 * this function approximates the fist order derivative 
 *
 * @param Y1 vector containing the n values of the image
 * @param Y2 vector containing the n values of the image offset by one value
 * @param X1 vector containing the n points where the image is valued
 * @param X2 vector containing the n points where the image is valued offset by one value
 *
 * @return derv vector containing the values of the derived function
 */
arma::mat Miscellaneous::deriv(arma::mat Y1, arma::mat Y2, arma::mat X1, arma::mat X2)
{
    arma::mat derv = (Y2 - Y1) / (X2 - X1);
    return derv;
}


void Miscellaneous::writeToTxt(arma::mat Z, std::string fileName)
{
    std::ofstream out("../out/" + fileName);
    if (out.is_open())
    {
        out << Z << std::endl;
        out.close();
    }
    std::cout << "The matrix has been written in " << fileName << std::endl;

}

void Miscellaneous::writeDataSet(arma::vec R, arma::vec Z, arma::mat rho)
{

    std::ofstream out("../out/result.csv");
    if (out.is_open())
    {
        std::string line = "";
        line = "Ri,Zj,Rhoij";
        out << line << std::endl;
        for (int i = 0; i < R.n_elem; i++)
        {
            for (int j = 0; j < Z.n_elem; j++)
            {
                line = to_string(R[i]) + "," + to_string(Z[j]) + "," + to_string(rho.at(i,j));
                out << line << std::endl;
            }
        }

    }
    std::cout << "Result dataset has been written in result.csv" << std::endl;
}

std::string Miscellaneous::cubeToDf3(const arma::cube &m)
{
  std::stringstream ss(std::stringstream::out | std::stringstream::binary);
  int nx = m.n_rows;
  int ny = m.n_cols;
  int nz = m.n_slices;
  ss.put(nx >> 8);
  ss.put(nx & 0xff);
  ss.put(ny >> 8);
  ss.put(ny & 0xff);
  ss.put(nz >> 8);
  ss.put(nz & 0xff);
  double theMin = 0.0;
  double theMax = m.max();
  for (uint k = 0; k < m.n_slices; k++)
  {
    for (uint j = 0; j < m.n_cols; j++)
    {
      for (uint i = 0; i < m.n_rows; i++)
      {
        uint v = 255 * (fabs(m(i, j, k)) - theMin) / (theMax - theMin);
        ss.put(v);
      }
    }
  }
  return ss.str();
}