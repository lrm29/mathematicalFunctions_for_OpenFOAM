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
    Foam::diagnostics

SourceFiles
    diagnostics.C

\*---------------------------------------------------------------------------*/

#ifndef diagnostics_H
#define diagnostics_H

#include "volFields.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//! \ingroup diagnostics
//! \brief Class to print diagnostics to screen.
//! \todo Maybe print some to file.
class diagnostics
{

	const fvMesh& mesh_;

public:


    // Constructors

        //- Construct from components
		diagnostics(const fvMesh& mesh);


    //- Destructor
    virtual ~diagnostics();


    // Member Functions

    //! Print minimum difference between two fields.
    void printMinDiffTwoFields(
    		const volScalarField& field1,
    		const volScalarField& field2
    ) const;

    //! Print mean, minimum and maximum values in a volScalarField.
    void meanMinMaxField(
    		const volScalarField& field
    ) const;
    //! Print mean, minimum and maximum values in a scalarField.
    void meanMinMaxField(
    		const scalarField& field,
    		const word& name
    ) const;

    //! Print mean, minimum and maximum values in a boundary.
    void meanMinMaxBoundary (
    		const volScalarField& field,
    		const word& boundary
    ) const;

    //! Catch negative values in a field.
    void catchNegativeValuesInField(
    		const volScalarField& field
    ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
