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
#include "EventAction.hh"

#include "SiPMHit.hh"
#include "SiliconPixelHit.hh"

#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4GenericMessenger.hh"
#include "g4root.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction() : G4UserEventAction() { DefineCommands(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event *) {
  fPrimariesPDG.clear();
  fPrimariesEnergy.clear();
  fPrimariesX.clear();
  fPrimariesY.clear();
  fPrimariesZ.clear();

  fSiHitsID.clear();
  fSiHitsX.clear();
  fSiHitsY.clear();
  fSiHitsZ.clear();
  fSiHitsEdep.clear();
  fSiHitsEdepNonIonising.clear();
  fSiHitsTOA.clear();
  fSiHitsTOAlast.clear();
  fSiHitsType.clear();

  fSiPMhitsID.clear();
  fSiPMhitsX.clear();
  fSiPMhitsY.clear();
  fSiPMhitsZ.clear();
  fSiPMhitsEdep.clear();
  fSiPMhitsEdepNonIonising.clear();
  fSiPMhitsTOA.clear();
  fSiPMhitsType.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event *event) {
  // sanity check
  if (event->GetNumberOfPrimaryVertex() == 0)
    return;

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleIColumn(0, event->GetEventID());
  // fill for all primary input particles
  for (G4int iVertex = 0; iVertex < event->GetNumberOfPrimaryVertex();
       iVertex++) {
    auto vertex = event->GetPrimaryVertex(iVertex);
    fPrimariesX.push_back(vertex->GetX0() / CLHEP::cm);
    fPrimariesY.push_back(vertex->GetY0() / CLHEP::cm);
    fPrimariesZ.push_back(vertex->GetZ0() / CLHEP::cm);
    for (G4int iParticle = 0; iParticle < vertex->GetNumberOfParticle();
         iParticle++) {
      auto particle = vertex->GetPrimary(iParticle);
      fPrimariesPDG.push_back(particle->GetPDGcode());
      fPrimariesEnergy.push_back(particle->GetTotalEnergy() / CLHEP::GeV);
    }
  }

  auto hce = event->GetHCofThisEvent();
  auto sdManager = G4SDManager::GetSDMpointer();
  G4int collId;

  // HGCAL EE + FH
  collId = sdManager->GetCollectionID("SiliconPixelHitCollection");
  auto hc = hce->GetHC(collId);
  if (!hc)
    return;
  double esumHGCAL = 0;
  double cogzHGCAL = 0;
  int NhitsHGCAL = 0;
  for (unsigned int i = 0; i < hc->GetSize(); ++i) {
    auto hit = static_cast<SiliconPixelHit *>(hc->GetHit(i));
    hit->Digitise(fHitTimeCut / CLHEP::ns, fToaThreshold / CLHEP::keV);

    if (hit->isValidHit()) {
      fSiHitsID.push_back(hit->ID());
      fSiHitsX.push_back(hit->GetX());
      fSiHitsY.push_back(hit->GetY());
      fSiHitsZ.push_back(hit->GetZ());
      fSiHitsEdep.push_back(hit->GetEdep());
      fSiHitsEdepNonIonising.push_back(hit->GetEdepNonIonizing());
      fSiHitsTOA.push_back(hit->GetTOA());
      fSiHitsTOAlast.push_back(hit->GetLastTOA());
      fSiHitsType.push_back(0);

      NhitsHGCAL++;
      esumHGCAL += hit->GetEdep() * CLHEP::keV / CLHEP::MeV;
      cogzHGCAL += hit->GetZ() * hit->GetEdep();
    }
  }
  if (esumHGCAL > 0)
    cogzHGCAL /= esumHGCAL;

  analysisManager->FillNtupleDColumn(23, esumHGCAL / CLHEP::GeV);
  analysisManager->FillNtupleDColumn(24, cogzHGCAL);
  analysisManager->FillNtupleIColumn(25, NhitsHGCAL);

  // AHCAL
  collId = sdManager->GetCollectionID("SiPMHitCollection");
  hc = hce->GetHC(collId);
  if (!hc)
    return;
  double esumAHCAL = 0;
  double cogzAHCAL = 0;
  int NhitsAHCAL = 0;
  for (unsigned int i = 0; i < hc->GetSize(); ++i) {
    auto hit = static_cast<SiPMHit *>(hc->GetHit(i));
    hit->Digitise(-1, 0);

    if (hit->isValidHit()) {
      fSiPMhitsID.push_back(hit->ID());
      fSiPMhitsX.push_back(hit->GetX());
      fSiPMhitsY.push_back(hit->GetY());
      fSiPMhitsZ.push_back(hit->GetZ());
      fSiPMhitsEdep.push_back(hit->GetEdep());
      fSiPMhitsEdepNonIonising.push_back(hit->GetEdepNonIonizing());
      fSiPMhitsTOA.push_back(hit->GetTOA());
      fSiPMhitsType.push_back(1);

      NhitsAHCAL++;
      esumAHCAL += hit->GetEdep() * CLHEP::keV / CLHEP::MeV;
      cogzAHCAL += hit->GetZ() * hit->GetEdep();
    }
  }
  if (esumAHCAL > 0)
    cogzAHCAL /= esumAHCAL;

  analysisManager->FillNtupleDColumn(26, esumAHCAL / CLHEP::GeV);
  analysisManager->FillNtupleDColumn(27, cogzAHCAL);
  analysisManager->FillNtupleIColumn(28, NhitsAHCAL);

  analysisManager->AddNtupleRow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::DefineCommands() {

  fMessenger = new G4GenericMessenger(this, "/HGCalTestbeam/hits/",
                                      "Primary generator control");

  // time cut command
  auto &timeCutCmd = fMessenger->DeclarePropertyWithUnit(
      "timeCut", "ns", fHitTimeCut,
      "Size of time window for hit digitalisation");
  timeCutCmd.SetParameterName("timeCut", true);
  timeCutCmd.SetRange("timeCut>=-1");
  timeCutCmd.SetDefaultValue("-1");

  // toa threshold command
  auto &toaThresholdCmd = fMessenger->DeclarePropertyWithUnit(
      "toaThreshold", "keV", fToaThreshold, "Threshold for TOA activation");
  toaThresholdCmd.SetParameterName("toaThreshold", true);
  toaThresholdCmd.SetRange("toaThreshold>=0");
  toaThresholdCmd.SetDefaultValue("0");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
