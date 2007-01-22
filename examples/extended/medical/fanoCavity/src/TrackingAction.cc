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
// $Id: TrackingAction.cc,v 1.2 2007-01-22 15:49:31 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "TrackingMessenger.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "HistoManager.hh"

#include "G4EmCalculator.hh"
#include "G4TrackingManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* det, RunAction* run,
                               SteppingAction* step, HistoManager* histo)
:detector(det),runAction(run),stepAction(step),histoManager(histo)
{
  matWall = 0;
  Zcav = 0.; 
  emCal = 0;
  first = true;
  killTrack = true;
  
  //create a messenger for this class
  trackMessenger = new TrackingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction()
{
  delete emCal;
  delete trackMessenger;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::SetKillTrack(G4bool flag)
{
  killTrack = flag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
 //get detector informations
 if (first) {
   matWall = detector->GetWallMaterial();
   Zcav    = 0.5*(detector->GetCavityThickness());
   emCal   = new G4EmCalculator();
   first   = false;
 }
  
 G4ParticleDefinition* particle = track->GetDefinition();
 G4bool charged = (particle->GetPDGCharge() != 0.);
 stepAction->TrackCharge(charged);
 G4int trackID = track->GetTrackID(); 
 if (trackID == 1 || !charged) return;
 
 //below, we have only charged secondaries
 //
 G4double energy = track->GetKineticEnergy(); 
 runAction->sumEsecond(energy);
 G4double position = (track->GetPosition()).z();
 // kill e- which cannot reach Cavity
 G4double safe = std::abs(position) - Zcav;
 G4double range = emCal->GetRangeFromRestricteDEDX(energy,particle,matWall);
 if (killTrack) {
   G4Track* aTrack = fpTrackingManager->GetTrack();   
   if (range < safe) aTrack->SetTrackStatus(fStopAndKill);
 }
   
 //histograms
 //
 histoManager->FillHisto(1,position);
 histoManager->FillHisto(2,energy);
 G4ThreeVector direction = track->GetMomentumDirection();
 histoManager->FillHisto(3,std::acos(direction.z()));      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track*)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

