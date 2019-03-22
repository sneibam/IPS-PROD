/**
 * @file Basis.cpp
 *
 * Fichier d'implémentation pour le calcul de la densité
 */

#include "Basis.h"


Basis::Basis()
{
}

Basis::Basis(double _br, double _bz, int N, double Q)
{
    bz = _bz;
    br = _br;

    int i = 0;
    int n, m;
    while (Basis::calc_nzMax(N,Q,i) >= 1)
    {
        i++;
    }
    mMax = i - 1;
    arma::ivec vm = arma::linspace<arma::ivec>(0, mMax-1, mMax);
    nMax = (mMax - vm - 1) / 2 + 1;

    n_zMax.zeros(vm.n_rows, (mMax - 1) / 2 + 1);

    for (m=0; m < mMax; m++)
        for (n=0; n < nMax[m]; n++)
        {
            n_zMax(m, n) = Basis::calc_nzMax(N, Q, m + 2*n + 1);
        }
}

double Basis::calc_nzMax(int N, double Q, int i)
{
    return (N + 2) * pow(Q, 2.0/3.0) + 0.5 - i * Q;
}


arma::vec Basis::zPart(arma::vec Z, int nz)
{
    arma::vec zpart;
    Poly poly;
    double factPart;
    double constPart;
    arma::vec expPart;

    
    factPart = 1/(std::sqrt(pow(2,nz) * std::sqrt(M_PI) * (double)Miscellaneous::factorial(nz)));
    arma::vec zfactPart =
    {
        7.511255444649424828587030047762276930524e-1,
        5.311259660135984572385365242537567693773e-1,
        2.655629830067992286192682621268783846887e-1,
        1.08415633823009689789878853323749022452e-1,
        3.833071493144388678758035489700725876944e-2,
        1.212123635259875349071310579655202864027e-2,
        3.499099535541983943542762678345281363805e-3,
        9.351736874441379770202142628079336715494e-4,
        2.337934218610344942550535657019834178873e-4,
        5.510563799824823821067738395272340242989e-5,
        1.232199525075784976451417893101893685019e-5,
        2.627058214392003194960140445607747615225e-6,
        5.362460124872919207572643465049275268564e-7,
        1.051664954522700642243897278328369487056e-7,
        1.987459951592227331357006655358164809117e-8,
        3.628588825417294648844562498457222314549e-9,
        6.414499411475746158307919785557463603948e-10,
        1.100077573465546081464447454369726996765e-10,
        1.833462622442576802440745757282878327941e-11,
        2.97426912202769525933866781584586998244e-12,
        4.702732399958400033923267262424369697562e-13,
    };

    //std::cout << "factPart done:" << std::endl << factPart << std::endl;
    constPart = (1 / std::sqrt(bz)) * zfactPart(nz);
    // constPart = (1 / std::sqrt(bz)) * factPart;
    //std::cout << "constPart done:" << std::endl << constPart << std::endl;
    expPart = (arma::exp(-(Z % Z) / (2 * bz * bz)));    
    //std::cout << "expPart done:" << std::endl << expPart << std::endl; 

    //std::cout << "bz: " << bz << endl;
    //std::cout << "Z: " << Z << endl;
    //std::cout << "nz: " << nz << endl;
    //std::cout << "Z/bz: " << Z/bz << endl;
    poly.calcHermite(nz, Z/bz);/*
    arma::vec hermitePart = poly.hermite(nz);
    hermitePart.print();
    std::cout << "Hermite done" << std::endl;*/

    // return ((1 / std::sqrt(bz)) * factPart) 
    //         % (arma::exp(-(Z % Z) / (2 * bz * bz)))
    //         % poly.hermite(nz);

    return constPart * expPart % poly.hermite(nz);
}



arma::vec Basis::rPart(arma::vec R, int m, int n)
{
    arma::vec rpart;
    Poly poly;

    poly.calcLaguerre(m+1, n+1, R % R / (br * br));

    return (sqrt((double) Miscellaneous::factorial(n)/ (double) Miscellaneous::factorial(n+m)) / (br * sqrt(M_PI))) 
        * (arma::exp(-(R % R) / (2 * br * br))) % pow(R / br, m) 
        % poly.laguerre(m,n);
}

arma::mat Basis::basisFunc(int m, int n, int nz, arma::vec Z, arma::vec R)
{
    return rPart(R, m, n) *  zPart(Z, nz).t();
}