//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
// $Id: HadrontherapyPhantomSD.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "HadrontherapyMatrix.hh"
#include "HadrontherapyAnalysisManager.hh"
#include "globals.hh"
#include <fstream>

HadrontherapyMatrix::HadrontherapyMatrix()
{  
  // Number of the voxels of the phantom
  numberVoxelX = 200;
  numberVoxelY = 200;
  numberVoxelZ = 200; 
 
  // Create the matrix
  matrix = new G4double[numberVoxelX*numberVoxelY*numberVoxelZ];
}

HadrontherapyMatrix::~HadrontherapyMatrix()
{
  delete[] matrix;
}
void HadrontherapyMatrix::Initialize()
{ 
  // Initialise the elemnts of the matrix to zero

  for(G4int i = 0; i < numberVoxelX; i++)
    {
      for(G4int j = 0; j < numberVoxelY; j++)
	{
	  for(G4int k = 0; k < numberVoxelZ; k++)

	    matrix[(i*numberVoxelY+j)*numberVoxelZ+k] = 0.;
	}
    }
}

void HadrontherapyMatrix::Fill(G4int i, G4int j, G4int k, 
			       G4double energyDeposit)
{
  if (matrix)
    matrix[(i*numberVoxelY+j)*numberVoxelZ+k] += energyDeposit;
  
  // Store the energy deposit in the matrix elemnt corresponding 
  // to the phantom voxel  
}

void HadrontherapyMatrix::TotalEnergyDeposit()
{
  // Store the information of the matrix in a ntuple and in 
  // a 1D Histogram

  G4int k;
  G4int j;
  G4int i;
  
  if (matrix)
    {  
	for(G4int l = 0; l < numberVoxelZ; l++) 
	  {
	    k = l;
	    
	    for(G4int m = 0; m < numberVoxelY; m++) 
	      { 
		j = m * numberVoxelZ + k; 
		
		for(G4int n = 0; n <  numberVoxelX; n++)
		  {
		    i =  n* numberVoxelZ * numberVoxelY + j;
		    if(matrix[i] != 0)
		      {	
					
#ifdef G4ANALYSIS_USE 	
			HadrontherapyAnalysisManager* analysis = 
			  HadrontherapyAnalysisManager::getInstance();
			analysis -> FillEnergyDeposit(n, m, k, matrix[i]);
			analysis -> BraggPeak(n, matrix[i]);
#endif
		      }
		  }       
	      }
	  }
    }
}
