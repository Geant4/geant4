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
// $Id: HadrontherapyMatrix.hh; May 2005
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

#ifndef HadrontherapyMatrix_H
#define HadrontherapyMatrix_H 1

#include "globals.hh"

// The information: energy deposit and position in the phantom
// is stored in a matrix

class HadrontherapyMatrix 
{
public:
  HadrontherapyMatrix();
  ~HadrontherapyMatrix();
 
  void Initialize(); 
  // All the elements of the matrix are initialised to zero
 
  void Fill(G4int i, G4int j, G4int k, G4double energyDeposit);
  // The matrix is filled with the energy deposit 
  // in the element corresponding to the voxel of the phantom where
  // the energy deposit was registered
 
  void TotalEnergyDeposit();
  // Store the information of the matrix in a ntuple and in 
  // a 1D Histogram

private:
  G4int numberVoxelX;
  G4int numberVoxelY;
  G4int numberVoxelZ;
  G4double* matrix;
};
#endif
