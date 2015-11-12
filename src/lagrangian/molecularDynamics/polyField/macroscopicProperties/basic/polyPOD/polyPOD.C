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
    polyPOD

Description

	Class to provide POD data smoothing functionality for the mdPODSmoothing
	post-processing application

\*----------------------------------------------------------------------------*/

#include "polyPOD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from IOdictionary and mesh
polyPOD::polyPOD
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict)
{
    // create new directory 

    fileName fieldPath(mesh_.time().path()/"PODsmoothed");
    fieldPath_ = fieldPath;
    
    if(isDir(fieldPath_) )
    {
        rmDir(fieldPath_);
    }

    mkDir(fieldPath_);
    
    const word fieldName = dict.lookup("fileName");
    
    fieldName_ = fieldName;
    
    nTimeSteps_ = readLabel(dict.lookup("nTimeSteps"));
    nBins_ = readLabel(dict.lookup("nBins"));
    
    massFlow_ = scalarRectangularMatrix(nTimeSteps_,nBins_);
    massFlowSmoother_ = scalarRectangularMatrix(nTimeSteps_,nBins_);
    
    for (label i=0; i<nTimeSteps_; i++)
    {
        for (label j=0; j<nBins_; j++)
        {
            massFlow_[i][j] = 0.0;
            massFlowSmoother_[i][j] = 0.0;
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyPOD::~polyPOD()
{}

void polyPOD::read()
{
    {
        Info << "Reading in matrix" << endl;
        
        word name = fieldName_;
        
        ifstream inFile;
        
        inFile.open(name.c_str());
        
        if (!inFile)
        {
            FatalErrorIn("pod::pod()")
            << "Unable to open file: " << name
            << exit(FatalError);
        }
        
        for (label i=0; i<nTimeSteps_; i++)
        {
            for (label j=0; j<nBins_; j++)
            {
                inFile >> massFlow_[i][j];
            }
        }
    }
    
    computePOD
    (
        massFlow_,
        nTimeSteps_,
        nBins_
    );
}
    
    

void polyPOD::write()
{

    Info << nl << "Writing smoothed matrix" << endl;
    
    OFstream file(fieldPath_/"smoothed.xy");
    
    if(file.good())
    {
        
        for (label i=0; i<nTimeSteps_; i++)
        {
            for (label j=0; j<nBins_; j++)
            {
                file << massFlowSmoother_[i][j] << " ";
            }
            
            file << endl;
        }        
    }
    else
    {
        FatalErrorIn("void write::write()")
        << "Cannot open file " << file.name()
        << abort(FatalError);
    }
}
    
void polyPOD::computePOD
(
    const scalarRectangularMatrix& data,
    const label& noSnapShots, //N
    const label& nBins // 200
)
{
    
    Info << "POD" << endl;
    
    scalarRectangularMatrix Q(noSnapShots, nBins);
    scalarRectangularMatrix Q2(noSnapShots, nBins);
    scalarField Qmean(nBins,0.0);
    
    for (label k=0; k<noSnapShots; k++)
    {
        for (label i=0; i<nBins; i++)
        {
            Q[k][i]=0.0;
            Q2[k][i]=0.0;
        }
    }
    
    Info<< "Size of the matrix:" << data.size() << endl;
    
    SVD svd(data, SMALL);
    
    //Sort singularValues in descending order and singularVectors accordingly
    
    scalarDiagonalMatrix S(svd.S().size());
    scalarRectangularMatrix U(svd.U().n(),svd.U().m());
    scalarRectangularMatrix V(svd.V().n(),svd.V().m());
    labelList indices(svd.S().size());
    
    SortableList<scalar> S_sorted;
    
    S_sorted=svd.S();
    S_sorted.reverseSort();
    indices=S_sorted.indices();
    reverse(indices);
    
    for (label k=0; k<svd.S().size(); k++)
    {
        S[k]=S_sorted[k];
    }
    
    for (label k=0; k<svd.U().n(); k++)
    {
        for (label i=0; i<svd.U().m(); i++)
        {
            U[k][i]=svd.U()[k][indices[i]];
            
        }
    }
    
    for (label k=0; k<svd.V().n(); k++)
    {
        for (label i=0; i<svd.V().m(); i++)
        {
            V[k][i]=svd.V()[k][indices[i]];
        }
    }
    
    multiply(Q, U, S, V.T());
    
    //Energy calculation.
    
    scalar sum_eigenVal=0.0;
    scalar energy=0.0;
    scalar energy_threshold=10.0;
    scalarList eigenVal(svd.S().size());
    eigenVal.setSize(svd.S().size(),0.0);
    
    for (label i=0; i<svd.S().size(); i++)
    {
        eigenVal[i]= pow(S[i],2);
        sum_eigenVal+= eigenVal[i];
    }
    
    Info<<"Sum of all eigenvalues:"<<sum_eigenVal<<endl;
    
    label modes=0;
    
    for (label i=0; i<svd.S().size(); i++)
    {
        energy=(eigenVal[i]/sum_eigenVal)*100;
        if (energy<energy_threshold) break;
        modes ++;
    }
    
    Info<<"Number of modes selected with energy test:"<<modes<<endl;
    
    scalarList svds(svd.S().size());
    svds.setSize(svd.S().size(),0.0);
    
    //POD approximation    
    scalarRectangularMatrix U_m(U.n(),modes);
    scalarRectangularMatrix V_m(U.m(),modes);
    scalarDiagonalMatrix S_m(modes);
    
    Info<<"*** Rebuilding the matrix ***"<<endl;
    Info<<"*** with "<<modes<<" mode(s) ***"<<endl;
    
    for (label m=0; m< modes; m++)
    {
        S_m[m]=S[m];
    }
    
    for (label j = 0; j < data.n(); j++)
    {
        for (label m=0; m< modes; m++)
        {
            U_m[j][m]=U[j][m];
        }
    }
    
    for (label i = 0; i < data.m(); i++)
    {
        for (label m=0; m< modes; m++)
        {
            V_m[i][m]=V[i][m];
        }
    }
    
    multiply(Q2, U_m, S_m, V_m.T());
    
    Info << "Setting smoother matrix" << endl;
    
    // set smoother matrix
    for(label i = 0; i < noSnapShots; ++i)
    {
        for(label j = 0; j < nBins; ++j)
        {
            massFlowSmoother_[i][j] = Q2[i][j];
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
