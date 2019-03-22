/**
 * @file Density.cpp
 *
 * Fichier pour le calcul de la densité
 */

#include "Density.h"


Density::Density()
{
    basis = Basis(1.935801664793151,      2.829683956491218,     14,     1.3);
}

/*
Density::Density(Basis _basis)
{
    basis = _basis;
}*/


/**
 *
 * This function computes the Density
 *
 * @param rVals values on axis R
 * @param zVals values on axis Z
 *
 * @return The density matrix
 */
arma::mat Density::calcDensity0(arma::vec rVals, arma::vec zVals)
{

    arma::mat rho;
    rho.load("rho.arma", arma::arma_ascii);

    int ia = 0;
    int ib = 0;
    arma::mat result = arma::zeros(rVals.n_rows, zVals.n_rows); // number of points on r- and z- axes
    for (int m = 0; m < basis.mMax; m++)
    {
        for (int n = 0; n < basis.nMax(m); n++)
        {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
            {
                for (int mp = 0; mp < basis.mMax; mp++)
                {
                    for (int np = 0; np < basis.nMax(mp); np++)
                    {
                        for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++)
                        {
                            arma::mat funcA = basis.basisFunc(m,  n,  n_z, zVals, rVals);
                            arma::mat funcB = basis.basisFunc(mp, np, n_zp, zVals, rVals);
                            result += funcA % funcB * rho(ia, ib); // mat += mat % mat * double
                            // std::cout << "Basis vector " << ia << " " << ib << ": m=" << mp << " n=" << np << " n_z=" << n_zp << std::endl;
                            ib++;
                        }
                    }
                }
                ib = 0;
                ia++;
            }            
        }
    }
    return result;
}


arma::mat Density::calcDensityOptimizedv1(arma::vec rVals, arma::vec zVals)
{

    arma::mat rho;
    rho.load("rho.arma", arma::arma_ascii);

    int ia = 0;
    int ib = 0;
    arma::mat result = arma::zeros(rVals.n_rows, zVals.n_rows); // number of points on r- and z- axes
    for (int m = 0; m < basis.mMax; m++)
    {
        for (int n = 0; n < basis.nMax(m); n++)
        {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
            {
                arma::mat funcA = basis.basisFunc(m,  n,  n_z, zVals, rVals);
                for (int mp = 0; mp < basis.mMax; mp++)
                {
                    for (int np = 0; np < basis.nMax(mp); np++)
                    {
                        for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++)
                        {
                            arma::mat funcB = basis.basisFunc(mp, np, n_zp, zVals, rVals);
                            result += funcA % funcB * rho(ia, ib); // mat += mat % mat * double
                            // std::cout << "Basis vector " << ia << " " << ib << ": m=" << mp << " n=" << np << " n_z=" << n_zp << std::endl;
                            ib++;
                        }
                    }
                }
                ib = 0;
                ia++;
            }            
        }
    }
    return result;
}


arma::mat Density::calcDensity1(arma::vec rVals, arma::vec zVals)
{

    arma::mat rho;
    rho.load("rho.arma", arma::arma_ascii);

    int ia = 0;
    int ib = 0;
    arma::mat result = arma::zeros(rVals.n_rows, zVals.n_rows); // number of points on r- and z- axes
    arma::imat sumN_zMax = arma::sum(basis.n_zMax, 1);

    for (int m = 0; m < basis.mMax; m++)
    {
        for (int n = 0; n < basis.nMax(m); n++)
        {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
            {
                for (int mp = 0; mp < basis.mMax; mp++)
                {
                    if (m != mp)
                    {
                        ib += sumN_zMax(mp);
                        continue;
                    }
                    for (int np = 0; np < basis.nMax(mp); np++)
                    {
                        for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++)
                        {
                            arma::mat funcA = basis.basisFunc(m,  n,  n_z, zVals, rVals);
                            arma::mat funcB = basis.basisFunc(mp, np, n_zp, zVals, rVals);
                            result += funcA % funcB * rho(ia, ib); // mat += mat % mat * double
                            // std::cout << "Basis vector " << ia << " " << ib << ": m=" << mp << " n=" << np << " n_z=" << n_zp << std::endl;
                            ib++;
                        }
                    }
                }
                ib = 0;
                ia++;
            }            
        }
    }
    return result;
}

arma::mat Density::calcDensity2(arma::vec rVals, arma::vec zVals)
{

    arma::mat rho;
    rho.load("rho.arma", arma::arma_ascii);

    int ia = 0;
    int ib = 0;
    arma::mat result = arma::zeros(rVals.n_rows, zVals.n_rows); // number of points on r- and z- axes
    arma::imat sumN_zMax = arma::sum(basis.n_zMax, 1);
    
    for (int m = 0; m < basis.mMax; m++)
    {
        for (int n = 0; n < basis.nMax(m); n++)
        {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
            {
                arma::mat funcA = basis.basisFunc(m,  n,  n_z, zVals, rVals);
                for (int mp = 0; mp < basis.mMax; mp++)
                {
                    if (m != mp)
                    {
                        ib += sumN_zMax(mp);
                        continue;
                    }
                    for (int np = 0; np < basis.nMax(mp); np++)
                    {
                        for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++)
                        {                            
                            arma::mat funcB = basis.basisFunc(mp, np, n_zp, zVals, rVals);
                            result += funcA % funcB * rho(ia, ib); // mat += mat % mat * double
                            // std::cout << "Basis vector " << ia << " " << ib << ": m=" << mp << " n=" << np << " n_z=" << n_zp << std::endl;
                            ib++;
                        }
                    }
                }
                ib = 0;
                ia++;
            }            
        }
    }
    return result;
}

arma::mat Density::calcDensity3(arma::vec rVals, arma::vec zVals)
{

    arma::mat rho;
    rho.load("rho.arma", arma::arma_ascii);

    int ia = 0;
    int ib = 0;
    arma::mat result = arma::zeros(rVals.n_rows, zVals.n_rows); // number of points on r- and z- axes
    arma::imat sumN_zMax = arma::sum(basis.n_zMax, 1);
    
    for (int m = 0; m < basis.mMax; m++)
    {
        for (int n = 0; n < basis.nMax(m); n++)
        {
            arma::vec rPartA = basis.rPart(rVals, m, n);
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
            {
                arma::vec zPartA = basis.zPart(zVals, n_z);
                arma::mat funcA = rPartA * zPartA.t();
                for (int mp = 0; mp < basis.mMax; mp++)
                {
                    if (m != mp)
                    {
                        ib += sumN_zMax(mp);
                        continue;
                    }
                    for (int np = 0; np < basis.nMax(mp); np++)
                    {
                        arma::vec rPartB = basis.rPart(rVals, mp, np);
                        for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++)
                        {     
                            //ib++;                       
                            arma::vec zPartB = basis.zPart(zVals, n_zp);
                            arma::mat funcB = rPartB * zPartB.t();
                            result += funcA % funcB * rho(ia, ib); // mat += mat % mat * double
                            // std::cout << "Basis vector " << ia << " " << ib << ": m=" << mp << " n=" << np << " n_z=" << n_zp << std::endl;
                            ib++;
                        }
                    }
                }
                ib = 0;
                ia++;
            }            
        }
    }
    return result;
}


arma::mat Density::calcDensity4(arma::vec rVals, arma::vec zVals)
{

    arma::mat rho;
    rho.load("rho.arma", arma::arma_ascii);

    int ia = 0;
    int ib = 0;
    arma::mat result = arma::zeros(rVals.n_rows, zVals.n_rows); // number of points on r- and z- axes
    arma::imat sumN_zMax = arma::sum(basis.n_zMax, 1);
    arma::mat sum = arma::zeros(rVals.n_rows, zVals.n_rows);
    
    for (int m = 0; m < basis.mMax; m++)
    {
        for (int n = 0; n < basis.nMax(m); n++)
        {
            arma::vec rPartA = basis.rPart(rVals, m, n);
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
            {
                arma::vec zPartA = basis.zPart(zVals, n_z);
                arma::mat funcA = rPartA * zPartA.t();
                sum = arma::zeros(rVals.n_rows, zVals.n_rows);
                for (int mp = 0; mp < basis.mMax; mp++)
                {
                    if (m != mp)
                    {
                        ib += sumN_zMax(mp);
                        continue;
                    }
                    for (int np = 0; np < basis.nMax(mp); np++)
                    {
                        arma::vec rPartB = basis.rPart(rVals, mp, np);
                        for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++)
                        {                          
                            arma::vec zPartB = basis.zPart(zVals, n_zp);
                            arma::mat funcB = rPartB * zPartB.t();
                            if (n_z > n_zp)
                            {
                                sum += 2 * funcB * rho(ia,ib);
                            }
                            else if (n_z == n_zp)
                            {
                                sum += funcB * rho(ia,ib);
                            }
                            ib++;
                        }
                    }
                }
                result += funcA % sum;
                ib = 0;
                ia++;
            }            
        }
    }
    return result;
}


// Deletion of the mp loop because when ma != mb the rho(a,b) = 0
arma::mat Density::calcDensity5(arma::vec rVals, arma::vec zVals)
{

    arma::mat rho;
    rho.load("rho.arma", arma::arma_ascii);

    int ia = 0;
    int ib = 0;
    int old_ib=0;
    arma::mat result = arma::zeros(rVals.n_rows, zVals.n_rows); // number of points on r- and z- axes
    arma::imat sumN_zMax = arma::sum(basis.n_zMax, 1);
    arma::mat sum = arma::zeros(rVals.n_rows, zVals.n_rows);
    
    for (int m = 0; m < basis.mMax; m++)
    {
        for (int n = 0; n < basis.nMax(m); n++)
        {
            arma::vec rPartA = basis.rPart(rVals, m, n);
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
            {
                arma::vec zPartA = basis.zPart(zVals, n_z);
                arma::mat funcA = rPartA * zPartA.t();
                sum = arma::zeros(rVals.n_rows, zVals.n_rows);
                ib = old_ib;
                for (int np = 0; np < basis.nMax(m); np++)
                {
                    arma::vec rPartB = basis.rPart(rVals, m, np);
                    for (int n_zp = 0; n_zp < basis.n_zMax(m, np); n_zp++)
                    {                          
                        arma::vec zPartB = basis.zPart(zVals, n_zp);
                        arma::mat funcB = rPartB * zPartB.t();
                        if (n_z > n_zp)
                        {
                            sum += 2 * funcB * rho(ia,ib); // mat += mat % mat * double
                        }
                        else if (n_z == n_zp)
                        {
                            sum += funcB * rho(ia,ib); // mat += mat % mat * double
                        }
                        ib++;
                    }
                }
                result += funcA % sum;
                ia++;
            }            
        }
        old_ib=ia;
    }
    return result;
}

arma::cube Density::transformIntoCylindricCoordinates(arma::vec R, arma::vec Y,arma::vec Z, arma::mat rho_plane)
{

    arma::cube rho3d(R.n_elem, Y.n_elem, Z.n_elem);
    for (int i = 0; i < R.n_elem; i++)
    {
        for (int j = 0; j < Y.n_elem; j++)
        {
            double r = sqrt(R(i)*R(i) + Y(i)*Y(i));
            int ind = 0;
            for (int k = 1; k < R.n_elem; k++)
            {
                if (abs(r - R[k]) < abs(r - R[ind]))
                {
                    ind = k;
                }
            }
            for (int k = 0; k < Z.n_elem; k++)
            {
                rho3d(i,j,k) = rho_plane(ind, k);
            }
        }
    }

    return rho3d;

}