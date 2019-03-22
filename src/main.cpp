/*! \mainpage <center style="color: #4182ea;">1D Quantum Harmonic Oscillator</center>
 *
 * In this project, we computed the solutions of the <span style="color: #4182ea;">1D Quantum Harmonic Oscillator</span> and we will check some of their properties.<br />
 * <strong>Schrödinger Equation</strong> :
 * <div style="font-size: 30px;"><center>\f$\hat{H}_{(z)}\psi_n(z) = E_n\psi_n(z)\f$</center></div>
 * with the 1D-Hamiltonian and 1D-momentum operators defined as <br/>
 * <div><center>\f$\hat{H}_{(z)}\equiv \frac{\hat{p}_{(z)}^2}{2m} + \frac{1}{2}m\omega^2\hat{z}^2, \hspace{5mm}\hat{p}_{(z)}\equiv -i\hbar\frac{\partial}{\partial z}.\f$</center></div>
 *
 * <div>Here are the solutions of the Schrödinger equation that we have to compute</div>
 * <center>\f$\psi_n(z) = \frac{1}{\sqrt{2^n n!}}\left(\frac{m\omega}{\pi\hbar}\right)^{1/4}e^{-\frac{m\omega z^2}{2\hbar}}H_n\left(\sqrt{\frac{m\omega}{\hbar}} . z\right).\f$</center>
 * The Hermmite polynomial \f$ H_n \f$ will be calculated using the recurence relation below :
 * <div>
 *  <center>\f$H_0(z) = 1\f$</center>
 * <center>\f$H_1(z) = 2z\f$</center>
 * <center>\f$\forall n\ge 1, H_{n+1}(z) = 2zH_n(z)-2nH_{n-1}(z).\f$</center>
 * </div>
 *
 * So after we computed the solutions we plotted the results using R.
 * <img src= "Rplot.png" width="800" />
 * <h3> Orthonormalité </h3>
 * <div style="margin: 15px;">
 * In order  to test our results we will compute the orthonomality. <br/>
 * <center>\f$\forall (m,n), \int \psi^*_m(z)\psi_n(z) dz = \delta_{mn}.\f$</center><br/>
 * We will use the quadrature rule of <strong>Carl Friedrich Gauss</strong> that allows us to approximate the definite integral of a function<br />
 * <center>\f$\int_{-1}^1 f(x) dx \simeq \sum_{i=0}^{n-1}w_if(x_i)\f$</center>
 * There are many quadrature but in our case we will use the <strong>Gauss-Hermite</strong> quadrature which looks as follow : <br />
 * <center><table style="background-color: rgb(0, 136, 192); border-collapse: collapse;" width="50%" border="1">
 *   <tr>
 *     <th>a</th>
 *     <th>b</th>
 *     <th>\f$w(x)\f$</th>
 *      <th>Quadrature</th>
 *     <th>Associated polynomial</th>
 *   </tr>
 *   <tr>
 *      <td style="background-color:white;" align="center">\f$-\infty\f$</td>
 *      <td style="background-color:white;" align="center">\f$+\infty\f$</td>
 *      <td style="background-color:white;" align="center">\f$\(e^{-x^2}\)\f$</td>
 *      <td style="background-color:white;" align="center"><strong>Gauss-Hermite</strong></td>
 *      <td style="background-color:white;" align="center">Hermite</td>
 *    </tr>
 * </table></center>
 * </div>
 */
#include <iostream>
#include <armadillo>
#include <fstream>
#include "Poly.h"
#include "Basis.h"
#include "Density.h"
#include "Miscellaneous.h"
using namespace arma;
using namespace std;

void test()
{
	Poly poly;
	arma::vec zVals, calcVals, targetVals, zVals2;
	zVals = {-3.1, -2.3, -1.0, -0.3, 0.1, 4.3, 9.2, 13.7};
	zVals2 = {-10.1, -8.4, -1.0, 0.0, 0.1, 4.3, 9.2, 13.7};
	poly.calcHermite(15, zVals2/2.829683956491218); // compute Hermite polynomials for n in {0 ... 5}
	calcVals = poly.hermite(15); // n = 4
	calcVals.print();
	// targetVals = {  1.02835360e+03,  2.05825600e+02, -2.00000000e+01,  7.80960000e+00,
	//                 1.15216000e+01,  4.59456160e+03,  1.10572154e+05,  5.54643458e+05};
}

void test2()
{
	cout << "laguerre" << endl;
	Poly poly;
	vec zVals, targetVals, calcVals;
	zVals = {0.1, 0.3, 1.2, 1.8, 2.0, 2.5, 7.1, 11.1};
	poly.calcLaguerre(6, 4, zVals);
	calcVals  = poly.laguerre(4, 2); // m = 4, n = 2
	//targetVals = {  14.405,  13.245,  8.52 ,  5.82 ,  5.,  3.125,  -2.395,  10.005};
}

void testZPart()
{
	Basis basis(1.935801664793151,      2.829683956491218,     14,     1.3);
	arma::vec z = {-10.1, -8.4, -1.0, 0.0, 0.1, 4.3, 9.2, 13.7};
	arma::vec res00 = { 7.64546544834383e-04,
			                    5.44886272162148e-03,
			                    4.19492564268520e-01,
			                    4.46522724110539e-01,
			                    4.46243982300708e-01,
			                    1.40736821086932e-01,
			                    2.26186220733178e-03,
			                    3.62929640195959e-06};
	z.print();
	arma::vec z0 = basis.zPart(z, 15);
}

void testDensity()
{
	arma::vec zVals, rVals;
	arma::mat basF, dens;
	Basis basis(1.935801664793151,      2.829683956491218,     14,     1.3);

	arma::vec z = {-10.1, -8.4, -1.0, 0.0, 0.1, 4.3, 9.2, 13.7};
	zVals = basis.zPart(z, 15);

	arma::vec r = {3.1, 2.3, 1.0, 0.0, 0.1, 4.3, 9.2, 13.7};
	rVals = basis.rPart(r, 8, 2);

	//basF = basis.basisFunc(8, 2, 15, z, r);
	Density density;
	dens = density.calcDensity0(rVals, zVals);
	dens.print();
	
	/*arma::mat rho;
    rho.load("rho.arma", arma::arma_ascii);
    std::cout << rho.n_rows << std::endl << rho.n_cols << endl;*/
}

void mainTest()
{
	Basis basis(1.935801664793151, 2.829683956491218, 14, 1.3);
	Density density;
	arma::wall_clock timer;
	arma::mat R;

	std::cout<< "Computing basis function . . ." << std::endl;
	arma::vec rVals= arma::linspace(-10,10,100);
    arma::vec zVals = arma::linspace(-10,10,100);

    arma::mat psi = basis.basisFunc(0, 0, 0, zVals, rVals);
    Miscellaneous::writeToTxt(rVals, "rVals_Psi000.txt");
    Miscellaneous::writeToTxt(zVals, "zVals_Psi000.txt");
	Miscellaneous::writeToTxt(psi, "Psi000.txt");

	psi = basis.basisFunc(0, 0, 1, zVals, rVals);
    Miscellaneous::writeToTxt(rVals, "rVals_Psi001.txt");
    Miscellaneous::writeToTxt(zVals, "zVals_Psi001.txt");
	Miscellaneous::writeToTxt(psi, "Psi001.txt");

	psi = basis.basisFunc(0, 1, 1, zVals, rVals);
    Miscellaneous::writeToTxt(rVals, "rVals_Psi011.txt");
    Miscellaneous::writeToTxt(zVals, "zVals_Psi011.txt");
	Miscellaneous::writeToTxt(psi, "Psi011.txt");

	psi = basis.basisFunc(1, 0, 1, zVals, rVals);
    Miscellaneous::writeToTxt(rVals, "rVals_Psi101.txt");
    Miscellaneous::writeToTxt(zVals, "zVals_Psi101.txt");
	Miscellaneous::writeToTxt(psi, "Psi101.txt");
	double n_seconds;
	std::cout<< "(version 0) Computing Density . . ." << std::endl;
	timer.tic();
	R = density.calcDensity0(rVals, zVals);
	n_seconds = timer.toc();
	std::cout << "(version 0) Number of seconds : " << n_seconds << std::endl;
	std::cout << "Mat Size : " << R.n_rows << "x" << R.n_cols << std::endl;
    Miscellaneous::writeToTxt(rVals, "rVals.txt");
    Miscellaneous::writeToTxt(zVals, "zVals.txt");
    Miscellaneous::writeToTxt(R, "plot2dv0.txt");
    Miscellaneous::writeDataSet(rVals, zVals, R);

	std::cout<< "(version 1: using kronecker delta) Computing Density . . ." << std::endl;
	timer.tic();
	R = density.calcDensity1(rVals, zVals);
	n_seconds = timer.toc();
	std::cout << "(version 1) Number of seconds : " << n_seconds << std::endl;
    //Miscellaneous::writeToTxt(rVals, "rVals.txt");
    //Miscellaneous::writeToTxt(zVals, "zVals.txt");
    //Miscellaneous::writeToTxt(R, "plot2dv1.txt");

    std::cout<< "(version 2: computing funcA before) Computing Density . . ." << std::endl;
	timer.tic();
	R = density.calcDensity2(rVals, zVals);
	n_seconds = timer.toc();
	std::cout << "(version 2) Number of seconds : " << n_seconds << std::endl;


	std::cout<< "(version 3: separating rPart & zPart) Computing Density . . ." << std::endl;
	timer.tic();
	R = density.calcDensity3(rVals, zVals);
	n_seconds = timer.toc();
	std::cout << "(version 3) Number of seconds : " << n_seconds << std::endl;

	//R.print();
	std::cout<< "(version 4: symmetry) Computing Density . . ." << std::endl;
	timer.tic();
	R = density.calcDensity4(rVals, zVals);
	n_seconds = timer.toc();
	std::cout << "(version 4) Number of seconds : " << n_seconds << std::endl;

	std::cout<< "(version 5: Removing the mp loop) Computing Density . . ." << std::endl;
	timer.tic();
	arma::mat R1 = density.calcDensity5(rVals, zVals);
	n_seconds = timer.toc();
	std::cout << "(version 5) Number of seconds : " << n_seconds << std::endl;

	std::cout << "Norm : " << arma::norm(R1-R) << std::endl;

}

void density3d()
{
	Density density;
	arma::vec rVals= arma::linspace(-10,10,32);
    arma::vec zVals = arma::linspace(-20,20,64);
    arma::vec yVals= arma::linspace(-10,10,32);
	arma::mat rho_plane = density.calcDensity0(rVals, zVals);
	arma::cube cube = density.transformIntoCylindricCoordinates(rVals, yVals, zVals, rho_plane);
	std::ofstream file;
    file.open("../out/results.df3");
    file << Miscellaneous::cubeToDf3(cube);
    file.close();


}

int main()
{
	//orthonormalityTest();
	//HermiteTest();
	//MatrixTest();
	//WnTest();
	//energyTest();
	//test();
	//testZPart();
	//testDensity();
	mainTest();
	// density3d();
	return 0;
}


