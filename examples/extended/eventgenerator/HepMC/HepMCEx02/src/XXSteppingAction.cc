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

// ====================================================================
//
//   XXSteppingAction.cc
//   $Id: XXSteppingAction.cc,v 1.1 2002-04-29 20:45:05 asaim Exp $
//
// ====================================================================

#include "XXSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

////////////////////////////////////
XXSteppingAction::XXSteppingAction()
////////////////////////////////////
{
}

/////////////////////////////////////
XXSteppingAction::~XXSteppingAction()
/////////////////////////////////////
{
}

//////////////////////////////////////////////////////////////
void XXSteppingAction::UserSteppingAction(const G4Step* astep)
  //////////////////////////////////////////////////////////////
{
  G4Track* atrack= astep-> GetTrack();

  if(atrack-> GetTrackStatus() != fAlive) return;
  if(atrack-> GetParentID() == 0) return;

  G4ParticleDefinition* particleType= atrack-> GetDefinition();
  if((particleType == G4MuonPlus::MuonPlusDefinition())
     || (particleType == G4MuonMinus::MuonMinusDefinition())) return;

  G4StepPoint* prestep= astep-> GetPreStepPoint();
  G4VPhysicalVolume* pv= prestep-> GetPhysicalVolume();
  G4String pvname= pv-> GetName();
  if(pvname=="BARREL_CAL_PV" || pvname=="ENDCAP_CAL_PV" ) {
    atrack-> SetTrackStatus(fKillTrackAndSecondaries);
  }
}

