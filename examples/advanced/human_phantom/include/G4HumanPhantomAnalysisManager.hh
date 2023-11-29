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

#ifndef G4HumanPhantomAnalysisManager_HH
#define G4HumanPhantomAnalysisManager_HH

#include "globals.hh"
#include "G4AnalysisManager.hh"

/*
Author: Susanna Guatelli, University of Wollongong, Australia

The class G4HumanPhantomAnalysisManager creates and manages ntuples

This class was developed following the extended Geant4 example analysis/AnaEx01

1 ntuple is created.
*/

// Define the total number of columns in the ntuple
const G4int MaxNtCol = 2;

class G4HumanPhantomAnalysisManager
{

public:
  G4HumanPhantomAnalysisManager();
  ~G4HumanPhantomAnalysisManager() = default;
  

  void book();
  // Create the output ROOT file 
  // Create the ntuple and histograms

  void FillNtupleWithEnergyDeposition(G4int,G4double);
  // Method to fill the ntuple with the energy deposition, integrated over a run, 
  // in each organ identified with an integer

  void save();
 // This method if called at the end of the run to store the 
 // results in the ROOT file

private:
    G4bool fFactoryOn; 
    G4int         fNtColId[MaxNtCol];
};
#endif



