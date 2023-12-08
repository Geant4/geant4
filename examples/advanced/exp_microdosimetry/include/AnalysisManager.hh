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
// Authors: Susanna Guatelli and Francesco Romano
// susanna@uow.edu.au, francesco.romano@ct.infn.it
//
 

#ifndef ANALYSISMANAGER_HH
#define ANALYSISMANAGER_HH 

#include "globals.hh"
#include "G4AnalysisManager.hh"

#include "AnalysisMessenger.hh"

class AnalysisMessenger;

// Define the total number of columns in the ntuple
const G4int MaxNtCol = 9;

class AnalysisManager
{ 

public:
   AnalysisManager(AnalysisMessenger* messenger);
  ~AnalysisManager();
  
  void book(G4bool addExtraNt); // booking the ROOT file

  void SetPrimaryEnergy(G4double energy); // Store the energy of the primary particles
  
  void StoreEnergyDeposition(G4double edep, G4double path, G4int eid);
  // Fill the ntuple with energy deposition and path length per event
  
  void FillSecondaries(G4int AA, G4double charge, G4double energy); 
  // Information about secondary particles
  
  void StoreSecondStageEnergyDeposition(G4double edep, G4int eid);
  // Fill the ntuple with energy and event ID of the second stage (telescope detector only)

  void finish();
  // Close the ROOT file with all the results stored in nutples 

private:
  G4bool factoryOn; 
  G4int         fNtColId[MaxNtCol];
  
  AnalysisMessenger* messenger;
  G4bool usingRoot;
  G4bool extraNt;

};

#endif




