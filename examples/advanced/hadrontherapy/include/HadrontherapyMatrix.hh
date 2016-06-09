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
// HadrontherapyMatrix.hh;
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy

#ifndef HadrontherapyMatrix_H
#define HadrontherapyMatrix_H 1

#include "globals.hh"
#include <vector>

// The information: energy deposit and position in the phantom
// is stored in a matrix

class HadrontherapyMatrix 
{
private:
  HadrontherapyMatrix(G4int voxelX, G4int voxelY, G4int voxelZ); //<--- this is supposed to be a singleton

public:

  ~HadrontherapyMatrix();
// Get object instance only
  static HadrontherapyMatrix* getInstance();
// Make & Get instance
  static HadrontherapyMatrix* getInstance(G4int nX, G4int nY, G4int nZ);
 
  void flush();
 
  void Initialize(); 
  // All the elements of the matrix are initialised to zero
   
  void Fill(G4int i, G4int j, G4int k, G4double energyDeposit);
  // The matrix is filled with the energy deposit 
  // in the element corresponding to the voxel of the phantom where
  // the energy deposit was registered
 
  void TotalEnergyDeposit();
  // Store the information of the matrix in a ntuple and in 
  // a 1D Histogram
  
  inline G4int Index(G4int i, G4int j, G4int k){ return (i * numberVoxelY + j) * numberVoxelZ + k; } 
  // Get a unique index from a three dimensional voxel information

private:

  static HadrontherapyMatrix* instance;
  G4int numberVoxelX;
  G4int numberVoxelY;
  G4int numberVoxelZ;
  G4double* matrix;
};
#endif
