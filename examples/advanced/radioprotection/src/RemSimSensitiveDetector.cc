//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: RemSimSensitiveDetector.cc,v 1.11 2005/11/23 09:22:11 guatelli Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Code developed by: S.Guatelli, guatelli@ge.infn.it
//
#include "RemSimSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#ifdef G4ANALYSIS_USE
#include "RemSimAnalysisManager.hh"
#endif

RemSimSensitiveDetector::RemSimSensitiveDetector(G4String name)
  :G4VSensitiveDetector(name)
{
}

RemSimSensitiveDetector::~RemSimSensitiveDetector(){;}

void RemSimSensitiveDetector::Initialize(G4HCofThisEvent*)
{

}

G4bool RemSimSensitiveDetector::ProcessHits(G4Step* aStep, 
                                            G4TouchableHistory* ROhist)
{
  G4double edep = aStep -> GetTotalEnergyDeposit();
  if(edep==0.) return false;

#ifdef G4ANALYSIS_USE
  RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
  G4int i = ROhist -> GetReplicaNumber();
  // Energy deposit in the phantom
  analysis -> energyDepositStore(i,edep/MeV);
  G4double xx = aStep -> GetPreStepPoint() -> GetPosition().x();
  G4double yy = aStep -> GetPreStepPoint() -> GetPosition().y();

  // Project the hits of primary and secondary particles
  // in the phantom in the plane x, y
  analysis -> particleShape(xx/cm, yy/cm);
 
  // Project the energy deposit of primary and secondary particles 
  // in the phantom in the plane x,y
  analysis -> energyDepShape(xx/cm,yy/cm, edep/MeV);

  // Energy deposit of secondary particles in the phantom
  if(aStep -> GetTrack() -> GetTrackID()!= 1)
    analysis -> SecondaryEnergyDeposit(i,edep/MeV);
#endif
  return true;
}

void RemSimSensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
}
