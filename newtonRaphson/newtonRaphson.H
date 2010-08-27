/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::newtonRaphson

SourceFiles
    newtonRaphson.C

\*---------------------------------------------------------------------------*/

#ifndef newtonRaphson_H
#define newtonRaphson_H

#include "dimensionedTypes.H"
#include <vector>
#include "scalarMatrices.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{



/*! \ingroup maths
 * \brief Class to perform the Newton-Raphson method of root finding.
 *
 * In order to use the newton raphson declare the following structure with the equations written into operator().
 * \verbatim
struct functionsToBeZeroed {

 	functionsToBeZeroed(){}
 	std::vector<Foam::scalar> operator()(const std::vector<Foam::scalar> x) {
 	    std::vector<Foam::scalar> y(1);
 	    y[0] = x[0]*x[0] -2.0;
 	    return y;
 	}

}; \endverbatim

In the main part of the program put the following with the relevant initial guess.

\verbatim
Foam::newtonRaphson nR;
bool check;
std::vector<Foam::scalar> initialGuess(1);
functionsToBeZeroed f;

initialGuess[0]=6.0;

nR.newt<functionsToBeZeroed>(initialGuess,check,f); \endverbatim

Also don't forget to put the following at the top of the header file:

\verbatim
#include "newtonRaphson.H" \endverbatim

 */
class newtonRaphson
{

	template <class T>
	void lnsrch(std::vector<Foam::scalar> &xold,
			    const Foam::scalar fold,
			    std::vector<Foam::scalar>  &g,
			    std::vector<Foam::scalar>  &p,
			    std::vector<Foam::scalar>  &x,
			    Foam::scalar &f,
			    const Foam::scalar stpmax,
			    bool &check,
			    T &func)
	{
	     const Foam::scalar ALF=1.0e-4, TOLX=1e-30;
	     Foam::scalar a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
	     Foam::scalar rhs1,rhs2,slope=0.0,sum=0.0,temp,test,tmplam;
	     int i,n=xold.size();
	     check=false;
	     for (i=0;i<n;i++) sum += p[i]*p[i];
	     sum=std::sqrt(sum);
	     if (sum > stpmax)
	          for (i=0;i<n;i++)
	            p[i] *= stpmax/sum;
	    for (i=0;i<n;i++)
	        slope += g[i]*p[i];
	    if (slope >= 0.0) throw("Roundoff problem in lnsrch.");
	    test=0.0;
	    for (i=0;i<n;i++) {
	        temp=std::abs(p[i])/Foam::max(std::abs(xold[i]),1.0);
	        if (temp > test) test=temp;
	    }
	    alamin=TOLX/test;
	    alam=1.0;
	    for (;;) {
	        for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
	        f=func(x);
	        if (alam < alamin) {
	            for (i=0;i<n;i++) x[i]=xold[i];
	            check=true;
	            return;
	        } else if (f <= fold+ALF*alam*slope) return;
	        else {
	            if (alam == 1.0)
	                tmplam = -slope/(2.0*(f-fold-slope));
	            else {
	                rhs1=f-fold-alam*slope;
	                rhs2=f2-fold-alam2*slope;
	                a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	                b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	                if (a == 0.0) tmplam = -slope/(2.0*b);
	                else {
	                    disc=b*b-3.0*a*slope;
	                    if (disc < 0.0) tmplam=0.5*alam;
	                    else if (b <= 0.0) tmplam=(-b+std::sqrt(disc))/(3.0*a);
	                    else tmplam=-slope/(b+std::sqrt(disc));
	                }
	                if (tmplam>0.5*alam)
	                    tmplam=0.5*alam;
	            }
	        }
	        alam2=alam;
	        f2 = f;
	        alam=Foam::max(tmplam,0.1*alam);
	    }
	}



	template <class T>
	struct NRfdjac {

	    const Foam::scalar EPS;
	    T &func;
	    NRfdjac(T &funcc) : EPS(1.0e-8),func(funcc) {}

	    Foam::scalarSquareMatrix operator() (std::vector<Foam::scalar> x, std::vector<Foam::scalar> fvec) {
	          int n=x.size();
	          Foam::scalarSquareMatrix df(n,n);
	          std::vector<Foam::scalar> xh=x;
	          for (int j=0;j<n;j++) {
	        	   Foam::scalar temp=xh[j];
	        	   Foam::scalar h=EPS*std::abs(temp);
	               if (h == 0.0) h=EPS;
	               xh[j]=temp+h;
	               h=xh[j]-temp;
	               std::vector<Foam::scalar> f=func(xh);
	               xh[j]=temp;
	               for (int i=0;i<n;i++)
	                   df[i][j]=(f[i]-fvec[i])/h;
	          }
	          return df;
	    }

	};



	template <class T>
	struct NRfmin {

	    std::vector<Foam::scalar> fvec;
	    T &func;
	    int n;
	    NRfmin(T &funcc) : func(funcc){}

	    Foam::scalar operator() (std::vector<Foam::scalar> x) {
	          n=x.size();
	          Foam::scalar sum=0;
	          fvec=func(x);
	          for (int i=0;i<n;i++) sum += Foam::sqr(fvec[i]);
	          return 0.5*sum;
	    }

	};


public:


	template <class T>
	void newt(std::vector<Foam::scalar>& x, bool &check, T &vecfunc) {

	    const int MAXITS=200;
	    const Foam::scalar TOLF=1.0e-8,TOLMIN=1.0e-12,STPMX=100.0;
	    const Foam::scalar TOLX=1e-30;
	    int i,j,its,n=x.size();
	    Foam::scalar den,f,fold,stpmax,sum,temp,test;
	    std::vector<Foam::scalar> g(n),p(n),xold(n);
	    Foam::scalarSquareMatrix fjac(n,n);
	    Foam::Field<scalar> solution(n);

	    NRfmin<T> fmin(vecfunc);
	    NRfdjac<T> fdjac(vecfunc);
	    std::vector<Foam::scalar> &fvec=fmin.fvec;
	    f=fmin(x);
	    test=0.0;

	    for (i=0;i<n;i++)
	         if (std::abs(fvec[i]) > test) test=std::abs(fvec[i]);
	    if (test < 0.01*TOLF) {
	         check=false;
	         return;
	    }
	    sum=0.0;
	    for (i=0;i<n;i++) sum += Foam::sqr(x[i]);
	    stpmax=STPMX*Foam::max(std::sqrt(sum),Foam::scalar(n));
	    for (its=0;its<MAXITS;its++) {
	         fjac=fdjac(x,fvec);
	         for (i=0;i<n;i++) {
	               sum=0.0;
	               for (j=0;j<n;j++) sum += fjac[j][i]*fvec[j];
	               g[i]=sum;
	         }
	         for (i=0;i<n;i++) xold[i]=x[i];
	         fold=f;
	         for (i=0;i<n;i++) solution[i] = -fvec[i];

	 		 //LUDecompose(fjac, pivotIndices);
	         //LUdcmp alu(fjac);
	         //alu.solve(p,p);
	         Foam::LUsolve<scalar>(fjac,solution);
	         for(int k=0;k<n;++k){
	        	 p[k] = solution[k];
	         }
	         lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin);

	         test=0.0;
	         for (i=0;i<n;i++)
	               if (std::abs(fvec[i]) > test) test=std::abs(fvec[i]);
	         if (test < TOLF) {
	               check=false;
	               return;
	         }
	         if (check) {
	               test=0.0;
	               den=Foam::max(f,0.5*n);
	               for (i=0;i<n;i++) {
	                   temp=std::abs(g[i])*Foam::max(std::abs(x[i]),1.0)/den;
	                   if (temp > test) test=temp;
	               }
	            check=(test < TOLMIN);
	            return;
	        }
	        test=0.0;
	        for (i=0;i<n;i++) {
	            temp=(std::abs(x[i]-xold[i]))/Foam::max(std::abs(x[i]),1.0);
	            if (temp > test) test=temp;
	        }
	        if (test < TOLX)
	            return;
	    }
	    throw("MAXITS exceeded in newt");
	}




    // Constructors

        //- Construct from components
		newtonRaphson();



    //- Destructor
    virtual ~newtonRaphson();


    // Member Functions



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
