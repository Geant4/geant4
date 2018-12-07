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
/// \file optical/LXe/src/LXeHistoManager.cc
/// \brief Implementation of the LXeHistoManager class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "LXeHistoManager.hh"
#include "G4UnitsTable.hh"

//#include<vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeHistoManager::LXeHistoManager()
  : fFileName("lxe")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeHistoManager::~LXeHistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeHistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in LXeHistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);    // enable inactivation of histograms

  // Define histogram indices, titles
  std::vector<std::pair<G4String, G4String> > histograms =
    { std::pair<G4String, G4String>("0", "dummy"),
      std::pair<G4String, G4String>("1", "hits per event"),
      std::pair<G4String, G4String>("2", "hits per event above threshold"),
      std::pair<G4String, G4String>("3", "scintillation photons per event"),
      std::pair<G4String, G4String>("4", "Cerenkov photons per event"),
      std::pair<G4String, G4String>("5", "absorbed photons per event"),
      std::pair<G4String, G4String>
                      ("6", "photons absorbed at boundary per event"),
      std::pair<G4String, G4String>
                      ("7", "energy deposition in scintillator per event"),
     };

  // Default values (to be reset via /analysis/h1/set command)
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (auto histogram : histograms) {
    G4int ih = analysisManager->
      CreateH1("h" + histogram.first, histogram.second, nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }
}
