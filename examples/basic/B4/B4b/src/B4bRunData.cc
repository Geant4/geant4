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
/// \file B4bRunData.cc
/// \brief Implementation of the B4bRunData class

#include "B4bRunData.hh"
#include "B4Analysis.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4bRunData::B4bRunData() 
 : G4Run(),
   fVolumeNames{ { "Absorber", "Gap" } }
{
  for ( auto& edep : fEdep ) { 
    edep = 0.; 
  }
  for ( auto& trackLength : fTrackLength ) {
    trackLength = 0.; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4bRunData::~B4bRunData()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4bRunData::FillPerEvent()
{
  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // accumulate statistic
  // in the order od the histograms, ntuple columns declarations
  G4int counter = 0;
  for ( auto edep : fEdep ) {
    analysisManager->FillH1(counter, edep);
    analysisManager->FillNtupleDColumn(counter++, edep);
  }
  for ( auto trackLength : fTrackLength ) {
    analysisManager->FillH1(counter, trackLength);
    analysisManager->FillNtupleDColumn(counter++, trackLength);
  }  
  analysisManager->AddNtupleRow();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4bRunData::Reset()
{ 
  for ( auto& edep : fEdep ) { 
    edep = 0.; 
  }
  for ( auto& trackLength : fTrackLength ) {
    trackLength = 0.; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
