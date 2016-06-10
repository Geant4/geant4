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
// $Id: B4bRunData.cc 69223 2013-04-23 12:36:10Z gcosmo $
//
/// \file B4bRunData.cc
/// \brief Implementation of the B4bRunData class

#include "B4bRunData.hh"
#include "B4Analysis.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4bRunData::B4bRunData() : G4Run()
{
  fVolumeNames[0] = "Absorber";
  fVolumeNames[1] = "Gap";
 
  for ( G4int i=0; i<kDim; i++) {
    fEdep[i] = 0.;
    fTrackLength[i] = 0.; 
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4bRunData::~B4bRunData()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4bRunData::FillPerEvent()
{
  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  //accumulate statistic
  //

  for ( G4int i=0; i<kDim; i++) {
    // fill histograms
    analysisManager->FillH1(i+1, fEdep[i]);
    analysisManager->FillH1(kDim+i+1, fTrackLength[i]);

    // fill ntuple
    analysisManager->FillNtupleDColumn(i, fEdep[i]);
    analysisManager->FillNtupleDColumn(kDim+i, fTrackLength[i]);
  }  

  analysisManager->AddNtupleRow();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4bRunData::Reset()
{ 
  for ( G4int i=0; i<kDim; i++) {
    fEdep[i] = 0.;
    fTrackLength[i] = 0.; 
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
