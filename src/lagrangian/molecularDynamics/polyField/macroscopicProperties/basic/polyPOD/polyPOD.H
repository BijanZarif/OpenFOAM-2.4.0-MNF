/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    pod

Description
    applied to a transient full MD simulation to compute mass flow rates and 
    density 

SourceFiles
    podI.H
    pod.C
    podIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyPOD_H
#define polyPOD_H

#include "Time.H"
#include "IOdictionary.H"
#include "IOstreams.H"
//#include "polyMesh.H"
#include "mdPoly.H"
#include "IFstream.H"

#include <iostream>
#include <iomanip>
#include <fstream>

#include "scalarMatrices.H"
#include "SVD.H"
#include "SortableList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pod Declaration
\*---------------------------------------------------------------------------*/

class polyPOD
{

private:

    // Private data

        const polyMesh& mesh_;
        const dictionary& dict_;
        fileName fieldPath_;
        word fieldName_;
        label nTimeSteps_;
        label nBins_;
        scalarRectangularMatrix massFlow_;
        scalarRectangularMatrix massFlowSmoother_;
        
    // Private Member Functions
    
        void computePOD
		(
		    const scalarRectangularMatrix& data,
            const label& noSnapShots, //N
            const label& nBins // 200
		);
        
public:

    // Constructors

        //- Construct from properties
    	polyPOD
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

    // Destructor

        ~polyPOD();

    // Member Functions
    
        void read();
        void write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
