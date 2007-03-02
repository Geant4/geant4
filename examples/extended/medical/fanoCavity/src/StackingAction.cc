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
// $Id: StackingAction.cc,v 1.2 2007-03-02 11:08:41 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"
#include "StackingMessenger.hh"

#include "G4Track.hh"
#include "G4EmCalculator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction(DetectorConstruction* det, RunAction* run,
                               HistoManager* histo)
:detector(det),runAction(run),histoManager(histo)
{
  matWall = 0;
  Zcav = 0.; 
  emCal = 0;
  first = true;
  killTrack  = true;
  
  //create a messenger for this class  
  stackMessenger = new StackingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{
  delete emCal;
  delete stackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* track)
{
 //get detector informations
 if (first) {
   matWall = detector->GetWallMaterial();
   Zcav    = 0.5*(detector->GetCavityThickness());
   emCal   = new G4EmCalculator();
   first   = false;
 }
 
 G4ClassificationOfNewTrack status = fUrgent;  

 //keep primary particle or neutral
 //
 G4ParticleDefinition* particle = track->GetDefinition();
 G4bool neutral = (particle->GetPDGCharge() == 0.);
 if ((track->GetParentID() == 0) || neutral) return status;

 //energy spectrum of charged secondaries
 //
  G4double energy = track->GetKineticEnergy();  
  runAction->sumEsecond(energy);
  
 // kill e- which cannot reach Cavity
 //
 G4double position = (track->GetPosition()).z();
 G4double safe = std::abs(position) - Zcav;
 G4double range = emCal->GetRangeFromRestricteDEDX(energy,particle,matWall);
 if (killTrack) {
   if (range < 0.8*safe) status = fKill;
 }

 //histograms
 //
 histoManager->FillHisto(1,position);
 histoManager->FillHisto(2,energy);
 G4ThreeVector direction = track->GetMomentumDirection();
 histoManager->FillHisto(3,std::acos(direction.z()));      
    
 return status;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
