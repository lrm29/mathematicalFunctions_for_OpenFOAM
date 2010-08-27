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
    Foam::incompleteGammaFunction

SourceFiles
    incompleteGammaFunction.C

\*---------------------------------------------------------------------------*/

#ifndef incompleteGammaFunction_H
#define incompleteGammaFunction_H

#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{



/*! \ingroup maths
 * \brief Class to calculate the incomplete gamma function.
 *
 * This is taken from numerical recipes. gammQ gives the upper incomplete gamma
 * function and gammP the lower incomplete gamma function.
 *
 * NOTE: Implementation has been altered slightly so that gammQ and gammP are
 * NOT normalised, as they have been in numerical recipes.
 *
 */
class incompleteGammaFunction
{

	static const int aSwitch = 100;
	static const scalar eps = 2.22e-16;
	static const scalar fpMin = 2.23e-308/2.22e-16;
	scalar gln;

	static const int ngau = 18;

	scalar gcf(const scalar a, const scalar x);
	scalar gSer(const scalar a, const scalar x);
	scalar gammln(const scalar xx);
	scalar gammPapprox(scalar a, scalar x, scalar psig);

public:


    // Constructors

        //- Construct from components
        incompleteGammaFunction();



    //- Destructor
    virtual ~incompleteGammaFunction();


    // Member Functions

	scalar gammP(const scalar a, const scalar x);
	scalar gammQ(const scalar a, const scalar x);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
