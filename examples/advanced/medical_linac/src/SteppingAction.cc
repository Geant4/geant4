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
// Code developed by
// Silvia Pozzi (1), silvia.pozzi@iss.it
// Barbara Caccia (1), barbara.caccia@iss.it
// Carlo Mancini Terracciano (2), carlo.mancini.terracciano@roma1.infn.it
// (1) Istituto Superiore di Sanita' and INFN Roma, Italy
// (2) Univ. La Sapienza and INFN Roma, Italy

#include "DetectorConstruction.hh"
#include "SteppingAction.hh"

#include <G4RunManager.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  DetectorConstruction* detector = (DetectorConstruction*)
    G4RunManager::GetRunManager() -> GetUserDetectorConstruction();
  
  G4double fFilterZ = detector -> GetFFilterZ();
  G4double fFilterRadius = detector -> GetFFilterRadius();
  
  G4String particleName;
  
  G4VPhysicalVolume* volume = aStep -> GetPostStepPoint() -> GetTouchable() -> GetVolume();
  
  const G4ThreeVector momentumDir = aStep -> GetTrack()->  GetMomentumDirection();
  if (momentumDir.z() < 0)
    {
		aStep->GetTrack()->SetTrackStatus(fStopAndKill);
    }
  
  if ( volume != nullptr )
    {
      G4String volumeName = volume -> GetName();
      
      if (volumeName == "Jaw1XPV" || volumeName == "Jaw2XPV"
	  || volumeName == "Jaw1YPV" || volumeName == "Jaw2YPV" )
	{
	  aStep->GetTrack()->SetTrackStatus(fStopAndKill);
	}
    }
  
  if (aStep -> GetPreStepPoint() -> GetTouchable() -> GetVolume() != nullptr )
    {
      if (aStep -> GetPreStepPoint() -> GetTouchable() -> GetVolume() -> GetName() == "TargetPV" ||
	  aStep -> GetPreStepPoint() -> GetTouchable() -> GetVolume() -> GetName() == "PriCollUpperPV" ||
	  aStep -> GetPreStepPoint() -> GetTouchable() -> GetVolume() -> GetName() == "PriCollMiddlePV" ||
	  aStep -> GetPreStepPoint() -> GetTouchable() -> GetVolume() -> GetName() == "PriCollLowererPV" ||
	  aStep -> GetPreStepPoint() -> GetTouchable() -> GetVolume() -> GetName() ==  "ffPV")
	{
	  if (std::abs(momentumDir.getX()/momentumDir.getZ()) > fFilterRadius/fFilterZ ||
	      std::abs(momentumDir.getY()/momentumDir.getZ()) > fFilterRadius/fFilterZ )
	    {
	      aStep->GetTrack()->SetTrackStatus(fStopAndKill);
	    }
	}
    }
}
