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

#ifndef BrachyAnalysisManager_HH
#define BrachyAnalysisManager_HH

#include "globals.hh"
#include "g4root.hh"

/*
Author: Susanna Guatelli

The class BrachyAnalysisManager creates and manages histograms and ntuples

This class was developed following the extended Geant4 example analysis/AnaEx01

1 ntuple is created.
*/
// Define the total number of histograms

// Define the total number of columns in the ntuple
const G4int MaxHisto = 1;
const G4int MaxNtCol = 4;

class BrachyAnalysisManager
{

public:
  BrachyAnalysisManager();
  ~BrachyAnalysisManager();
  

  void book();
  // Create the output ROOT file 
  // Create the ntuple and histograms

  void FillNtupleWithEnergyDeposition(G4double,G4double,G4double,G4double);
  // Method to fill the ntuple with the energy deposition, integrated over a run, in each voxel
  // of the scoring mesh

  void FillPrimaryParticleHistogram(G4double);
  // Energy spectrum of primary particles

  void save();
 // This method if called at the end of the run to store the 
 // results in the ROOT file

private:
    G4bool factoryOn; 
    G4int         fHistId[MaxHisto];
    G4int         fNtColId[MaxNtCol];
    G4AnaH1*      primaryParticleSpectrum;
};
#endif



