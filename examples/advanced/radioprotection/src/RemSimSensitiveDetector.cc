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
// $Id: RemSimSensitiveDetector.cc,v 1.15 2006-06-29 16:24:19 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

  i = ROhist -> GetReplicaNumber();

#ifdef G4ANALYSIS_USE
  RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
  
  // Energy deposit in the phantom
  analysis -> energyDepositStore(i,edep/MeV);
  G4double xx = aStep -> GetPreStepPoint() -> GetPosition().x();
  G4double yy = aStep -> GetPreStepPoint() -> GetPosition().y();
  //  G4double zz = aStep -> GetPreStepPoint() -> GetPosition().z();

  // Project the hits of primary and secondary particles
  // in the phantom in the plane x, y
  G4String volume = aStep -> GetPreStepPoint()-> GetPhysicalVolume() -> GetName();
 
  if ( (xx < 150. *cm && xx > -150. * cm) && (yy < 150. *cm && yy > -150. * cm)  )
  { 
  analysis -> particleShape(xx/cm, yy/cm);
  // 
  // Project the energy deposit of primary and secondary particles 
  // in the phantom in the plane x,y
  analysis -> energyDepShape(xx/cm,yy/cm, edep/MeV);
  }

  // Energy deposit of secondary particles in the phantom
  if(aStep -> GetTrack() -> GetTrackID()!= 1)
    analysis -> SecondaryEnergyDeposit(i,edep/MeV);

			     
#endif
  return true;
}

void RemSimSensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
}
