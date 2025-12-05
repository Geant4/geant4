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
/// \file GB03EventAction.cc
/// \brief Implementation of the GB03EventAction class

#include "GB03EventAction.hh"

#include "GB03HodoHit.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4THitsMap.hh"
#include "G4UnitsTable.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03EventAction::BeginOfEventAction(const G4Event* /*event*/) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03EventAction::EndOfEventAction(const G4Event* event)
{
  // Get hist collections IDs

  auto* sdman = G4SDManager::GetSDMpointer();
  auto hodoCID = sdman->GetCollectionID("HodoSD/Leaving");

  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // fill histograms

  auto* hodoCollection =
    static_cast<GB03HodoHitsCollection*>(event->GetHCofThisEvent()->GetHC(hodoCID));

  for (std::size_t i = 0; i < hodoCollection->entries(); i++) {
    auto id = (*hodoCollection)[i]->GetId();
    auto w = (*hodoCollection)[i]->GetWeight();
    auto Ek = (*hodoCollection)[i]->GetEkin();
    auto xpos = (*hodoCollection)[i]->GetPos().x();
    if (id == 2112) {  // neutron
      analysisManager->FillH1(0, Ek, w);  // Ekin (Mev)
      analysisManager->FillH1(2, xpos / cm, w);  // n hits - x pos
    }
    else if (id == 22) {  // gamma
      analysisManager->FillH1(1, Ek, w);  // Ekin (Mev)
      analysisManager->FillH1(3, xpos / cm, w);  // gamma hits - x pos
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
