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
/// \file optical/wls/src/WLSEventAction.cc
/// \brief Implementation of the WLSEventAction class
//
//

#include "WLSEventAction.hh"

#include "WLSEventActionMessenger.hh"
#include "WLSPhotonDetHit.hh"
#include "WLSRun.hh"
#include "WLSRunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VVisManager.hh"
#include "Randomize.hh"

// Purpose: Accumulates statistics regarding hits
//          in the PhotonDet detector

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSEventAction::WLSEventAction()
  : fVerboseLevel(0)
{
  fMPPCCollID = 0;

  fEventMessenger = new WLSEventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSEventAction::~WLSEventAction() { delete fEventMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSEventAction::BeginOfEventAction(const G4Event*)
{
  fNTIR        = 0;
  fNExiting    = 0;
  fEscapedEnd  = 0;
  fEscapedMid  = 0;
  fBounce      = 0;
  fWLSBounce   = 0;
  fClad1Bounce = 0;
  fClad2Bounce = 0;
  fReflected   = 0;
  fEscaped     = 0;
  fMirror      = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSEventAction::EndOfEventAction(const G4Event* evt)
{
  // Get Hits from the detector if any
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String colName   = "PhotonDetHitCollection";
  fMPPCCollID        = SDman->GetCollectionID(colName);

  G4HCofThisEvent* HCE               = evt->GetHCofThisEvent();
  WLSPhotonDetHitsCollection* mppcHC = nullptr;

  // Get the hit collections
  if(HCE)
  {
    if(fMPPCCollID >= 0)
    {
      mppcHC = (WLSPhotonDetHitsCollection*) (HCE->GetHC(fMPPCCollID));
    }
  }

  // Get hit information about photons that reached the detector in this event
  G4int n_hit = 0;
  if(mppcHC)
  {
    n_hit = mppcHC->entries();
  }

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillH1(2, mppcHC->entries());
  for (size_t i = 0; i < mppcHC->entries(); ++i) {
    auto pdHit = (*mppcHC)[i];
    analysisManager->FillH1(0, pdHit->GetEnergy());
    analysisManager->FillH1(1, pdHit->GetArrivalTime());
  }

  if(fVerboseLevel > 1)
  {
    G4cout << "-------------------------------------" << G4endl
           << " In this event, number of:" << G4endl
           << "  TIR:           " << fNTIR << G4endl
           << "  Exiting:       " << fNExiting << G4endl
           << "  Escaped Mid:   " << fEscapedMid << G4endl
           << "  Escaped End:   " << fEscapedEnd << G4endl
           << "  Bounced:       " << fBounce << G4endl
           << "  WLS Bounce:    " << fWLSBounce << G4endl
           << "  Clad1 Bounce:  " << fClad1Bounce << G4endl
           << "  Clad2 Bounce:  " << fClad2Bounce << G4endl
           << "  Reflected:     " << fReflected << G4endl
           << "  Escaped:       " << fEscaped << G4endl
           << "  Mirror:        " << fMirror << G4endl
           << "  Detector hit:  " << n_hit << G4endl;
  }

  WLSRun* run = static_cast<WLSRun*>(
    G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  run->AddTIR(fNTIR);
  run->AddExiting(fNExiting);
  run->AddEscapedEnd(fEscapedEnd);
  run->AddEscapedMid(fEscapedMid);
  run->AddBounce(fBounce);
  run->AddClad1Bounce(fClad1Bounce);
  run->AddClad2Bounce(fClad2Bounce);
  run->AddReflected(fReflected);
  run->AddEscaped(fEscaped);
  run->AddMirror(fMirror);
  run->AddDetectorHits(n_hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int WLSEventAction::GetEventNo()
{
  return fpEventManager->GetConstCurrentEvent()->GetEventID();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSEventAction::SetEventVerbose(G4int level) { fVerboseLevel = level; }
