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
// $Id: G4Step.cc,v 1.11 2010-10-30 07:49:08 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
//  G4Step.cc
//
//  Description:
//    This class represents the Step of a particle tracked.
//    It includes information of 
//      1) List of Step points which compose the Step,
//      2) static information of particle which generated the 
//         Step, 
//      3) trackID and parent particle ID of the Step,
//      4) termination condition of the Step,
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
// ---------------------------------------------------------------

#include "G4Step.hh"

////////////////
G4Step::G4Step()
////////////////
  :  fTotalEnergyDeposit(0.0),fNonIonizingEnergyDeposit(0.0),
     fStepLength(0.), fpTrack(0), 
     fpSteppingControlFlag(NormalCondition),fSecondary(0),
     fpVectorOfAuxiliaryPointsPointer(0)
{
  fpPreStepPoint  = new G4StepPoint();
  fpPostStepPoint = new G4StepPoint();

  fFirstStepInVolume =false;
  fLastStepInVolume = false;
}

/////////////////
G4Step::~G4Step()
/////////////////
{
  delete fpPreStepPoint;
  delete fpPostStepPoint;
}


/////////////////
G4ThreeVector G4Step::GetDeltaMomentum() const
/////////////////
{ 
  static G4bool isFirstTime = true;
  if (isFirstTime) {
    isFirstTime = false;
#ifdef G4VERBOSE
    G4Exception( "G4Step::GetDeltaMomentum()","Warning", JustWarning, 
		 "This method is obsolete and will be removed soon");
#endif
  }

  return fpPostStepPoint->GetMomentum()
    - fpPreStepPoint->GetMomentum(); 
}

/////////////////
G4double G4Step::GetDeltaEnergy() const
  /////////////////
{ 
  static G4bool isFirstTime = true;
  if (isFirstTime) {
    isFirstTime = false;
#ifdef G4VERBOSE
    G4Exception( "G4Step::GetDeltaEnergy()","Warning", JustWarning, 
		 "This method is obsolete and will be removed soon");
#endif
  }

  return fpPostStepPoint->GetKineticEnergy()
    - fpPreStepPoint->GetKineticEnergy(); 
}
