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
    Foam::numericalIntegration

SourceFiles
    numericalIntegration

\*---------------------------------------------------------------------------*/

#ifndef numericalIntegration_H
#define numericalIntegration_H

#include "dimensionedTypes.H"
#include <vector>
#include "scalarMatrices.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{



/*! \ingroup maths
 * \brief Class to perform numerical integration.
 */
class numericalIntegration
{

	struct Quadrature{
	    int n;
	    virtual double next() = 0;
	};

	/*!
	 * Midpoint quadrature.
	 *
	 * This is an open formula making it suitable for calculating improper integrals.
	 */
	template <class T>
	struct Midpnt : Quadrature {

	    double a,b,s;
	    T &funk;
	    Midpnt(T &funcc, const double aa, const double bb) :
	         funk(funcc), a(aa), b(bb) {n=0;}
	    double next(){
	         int it,j;
	         double x,tnm,sum,del,ddel;
	         n++;
	         if (n == 1) {
	              return (s=(b-a)*func(0.5*(a+b)));
	         } else {
	              for(it=1,j=1;j<n-1;j++) it *= 3;
	              tnm=it;
	              del=(b-a)/(3.0*tnm);
	              ddel=del+del;
	              x=a+0.5*del;
	              sum=0.0;
	              for (j=0;j<it;j++) {
	                   sum += func(x);
	                   x += ddel;
	                   sum += func(x);
	                   x += del;
	              }
	              s=(s+(b-a)*sum/tnm)/3.0;
	              return s;
	         }
	    }
	    virtual double func(const double x) {return funk(x);}
	};

	/*!
	 * Trapezoidal quadrature.
	 */
	template<class T>
	struct Trapzd : Quadrature {

		double a,b,s;
	    T &func;
	    Trapzd() {};
	    Trapzd(T &funcc, const double aa, const double bb) :
	         func(funcc), a(aa), b(bb) {n=0;}
	    double next() {
	    	double x,tnm,sum,del;
	         int it,j;
	         n++;
	         if (n == 1) {
	              return (s=0.5*(b-a)*(func(a)+func(b)));
	         } else {
	              for (it=1,j=1;j<n-1;j++) it <<= 1;
	              tnm=it;
	              del=(b-a)/tnm;
	              x=a+0.5*del;
	              for (sum=0.0,j=0;j<it;j++,x+=del) sum += func(x);
	              s=0.5*(s+(b-a)*sum/tnm);
	              return s;
	         }
	    }

	};


public:

	/*!
	 * Driver routine for trapezoidal integration.
	 */
	template<class T>
	double qtrap(T &func, const double a, const double b, const double eps=1.0e-6) {

	     const int JMAX=20;

	     double s,olds=0.0;
	     Trapzd<T> t(func,a,b);
	     for (int j=0;j<JMAX;j++) {
	         s=t.next();
	         if (j > 5)
	              if (std::abs(s-olds) < eps*std::abs(olds) ||
	                   (s == 0.0 && olds == 0.0)) return s;
	         olds=s;
	     }
	     throw("Too many steps in routine qtrap");

	}

	/*!
	 * Driver routine for trapezoidal integration with fixed step size.
	 */
	template<class T>
	double qtrapfixed(T &func, const double a, const double b, const double m=5) {

		double s=0.0;

		Trapzd<T> t(func,a,b);
		for(int j=1;j<=m+1;j++) s=t.next();

		return s;

	}

	/*!
	 * Driver routine for trapezoidal integration.
	 */
	template<class T>
	double qmid(T &func, const double a, const double b, const double eps=1.0e-6) {

	     const int JMAX=20;

	     double s,olds=0.0;
	     Midpnt<T> t(func,a,b);
	     for (int j=0;j<JMAX;j++) {
	         s=t.next();
	         if (j > 5)
	              if (std::abs(s-olds) < eps*std::abs(olds) ||
	                   (s == 0.0 && olds == 0.0)) return s;
	         olds=s;
	     }
	     throw("Too many steps in routine qmid");

	}


    // Constructors

        //- Construct from components
		numericalIntegration();



    //- Destructor
    virtual ~numericalIntegration();


    // Member Functions



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
