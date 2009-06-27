//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: HadrontherapyMatrix.cc;
// Last modified: G.A.P.Cirrone
//
// This file create a file where information on dose deposited in each matrix voxel 
// is stored
//
// See more at: http://workgroup.lngs.infn.it/geant4lns
//
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
HadrontherapyMatrix* HadrontherapyMatrix::instance = 0;
HadrontherapyMatrix* HadrontherapyMatrix::getInstance()
{
  if (instance == 0) instance = new HadrontherapyMatrix;
  return instance;
}
HadrontherapyMatrix::HadrontherapyMatrix()
{  
// Number of the voxels of the phantom
// For Y = Z = 1 the phantom is dived in slices (and not in voxels)
// orthogonal to the beam axis
  numberVoxelX = 200;
  numberVoxelY = 1;
  numberVoxelZ = 1; 
 
  // Create the matrix
  matrix = new G4double[numberVoxelX*numberVoxelY*numberVoxelZ];
}

HadrontherapyMatrix::~HadrontherapyMatrix()
{
  delete[] matrix;
}

void HadrontherapyMatrix::flush(){
	if(matrix)
		for(int i=0;i<numberVoxelX*numberVoxelY*numberVoxelZ;i++){
		matrix[i] = 0;
		}
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
      std::ofstream ofs;    	
      ofs.open("DoseDistribution.out"); 
      
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
			ofs<< n <<'\t'<< m <<'\t'<<
			  k<<'\t'<<matrix[i]<<G4endl;
		       	
#ifdef ANALYSIS_USE
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
