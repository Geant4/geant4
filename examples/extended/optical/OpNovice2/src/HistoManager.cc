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
/// \file optical/OpNovice2/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("opnovice2")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);    // enable inactivation of histograms

  // Define histogram indices, titles
  G4int maxHisto = 15;
  G4String id[] = { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                   "10","11","12","13","14","15","16","17","18","19" };
 
  // TODO change throughout code
  G4String title[] = {
      "dummy",                                        // 0
      "Cerenkov spectrum",                            // 1
      "scintillation spectrum",                       // 2
      "scintillation photons creation time",          // 3
      "WLS absorption spectrum",                      // 4
      "WLS emission spectrum",                        // 5
      "WLS emission time",                            // 6
      "WLS2 absorption spectrum",                     // 7
      "WLS2 emission spectrum",                       // 8
      "WLS2 emission time",                           // 9
      "boundary process status",                      //10
      "X momentum dir of backward-going photons",     //11
      "Y momentum dir of backward-going photons",     //12
      "Z momentum dir of backward-going photons",     //13
      "X momentum dir of forward-going photons",      //14
      "Y momentum dir of forward-going photons",      //15
      "Z momentum dir of forward-going photons",      //16
      "X momentum dir of Fresnel-refracted photons",  //17
      "Y momentum dir of Fresnel-refracted photons",  //18
      "Z momentum dir of Fresnel-refracted photons",  //19
  };
 
  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  for (G4int k=0; k <= maxHisto; ++k) {
      G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
      analysisManager->SetH1Activation(ih, false);
  }
}
