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

\*---------------------------------------------------------------------------*/

#include "diagnostics.H"

// * * * * * * * * * * * * * * * * Constructors* * * * * * * * * * * * * * * //

Foam::diagnostics::diagnostics(const fvMesh& mesh)
:
	mesh_(mesh)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diagnostics::~diagnostics()
{}

// * * * * * * * * * * * * * * * * Member Functions* * * * * * * * * * * * * //

/*!
 * Print to screen the minimum difference between two fields. Does not include boundaries.
 *
 * \param[in] volScalarField Field1
 * \param[in] volScalarField Field2
 */
void
Foam::diagnostics::printMinDiffTwoFields(const volScalarField& field1, const volScalarField& field2)
const
{
    Info << "Minimum difference between fields "
		 << field1.name()
		 << " and "
		 << field2.name()
		 << " = "
		 << Foam::min(Foam::mag(field1.internalField() - field2.internalField()))
    	 << endl;
}

/*!
 * Print to screen the mean, minimum and maximum of a volScalarField (includes boundaries).
 *
 * \param[in] volScalarField field
 */
void
Foam::diagnostics::meanMinMaxField(const volScalarField& field)
const
{
	Info << "Stats for field "
	     << field.name() << ":"
	     << " Weighted Mean = " << field.weightedAverage(mesh_.V()).value()
	     << " Min = " << Foam::gMin(field)
	     << " Max = " << Foam::gMax(field)
	     << endl;
}

/*!
 * Print to screen the mean, minimum and maximum of a scalarField.
 *
 * \param[in] scalarField field
 * \param[in] word	Name of the field (supplied by user)
 */
void
Foam::diagnostics::meanMinMaxField(const scalarField& field, const word& name)
const
{
	Info << "Stats for field "
	     << name << ":"
	     << " Mean = " << Foam::gAverage(field)
	     << " Min = " << Foam::gMin(field)
	     << " Max = " << Foam::gMax(field)
	     << endl;
}

/*!
 * Print to screen the mean, minimum and maximum of a boundary.
 *
 * \param[in] volScalarField field
 * \param[in] word boundary
 *
 * NOTE: Does not currently work in parallel.
 *
 */
void
Foam::diagnostics::meanMinMaxBoundary(const volScalarField& field, const word& boundary)
const
{
	if(!Pstream::parRun()){
		forAll(field.boundaryField(), patchI)
		{
			if(mesh_.boundaryMesh()[patchI].name() == boundary)
			{
				Info << "Stats for boundary "
					 << boundary << ":"
					 << " Mean = " << Foam::average(field.boundaryField()[patchI])
					 << " Min = " << Foam::min(field.boundaryField()[patchI])
					 << " Max = " << Foam::max(field.boundaryField()[patchI])
					 << endl;
			}
		}
	} else {
		Info << "Parallel run: meanMinMaxBoundary may not currently work in parallel" << endl;

		scalar boundaryFieldSize,boundaryFieldSum;

		forAll(field.boundaryField(), patchI)
		{
			if(mesh_.boundaryMesh()[patchI].name() == boundary
			   && field.boundaryField()[patchI].size() != 0)
			{

				boundaryFieldSize = field.boundaryField()[patchI].size();
				boundaryFieldSum = sum(field.boundaryField()[patchI]);
			}
		}

		reduce(boundaryFieldSize,sumOp<scalar>());
		reduce(boundaryFieldSum,sumOp<scalar>());


				Info << "Stats for boundary "
					 << boundary << ":"
					 << " Mean = " << boundaryFieldSum/boundaryFieldSize
					 //<< " Min = " << boundaryFieldSum/
						 //<< " Max = " << Foam::reduce(field.boundaryField()[patchI])
					 << endl;

	}
}

/*!
 * Catch negative values in a field.
 *
 * \param[in] volScalarField field
 */
void
Foam::diagnostics::catchNegativeValuesInField(const volScalarField& field)
const
{
	if(Foam::min(field).value() < 0.0)
	{
		Pout << "Caught: Value in field " << field.name() << " less than zero: "
			 << endl;
	}
}

// ************************************************************************* //
