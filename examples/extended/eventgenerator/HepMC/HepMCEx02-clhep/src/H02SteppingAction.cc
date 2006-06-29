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
// ====================================================================
//
//   H02SteppingAction.cc
//   $Id: H02SteppingAction.cc,v 1.3 2006-06-29 17:14:16 gunter Exp $
//
// ====================================================================

#include "H02SteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

////////////////////////////////////
H02SteppingAction::H02SteppingAction()
////////////////////////////////////
{
}

/////////////////////////////////////
H02SteppingAction::~H02SteppingAction()
/////////////////////////////////////
{
}

//////////////////////////////////////////////////////////////
void H02SteppingAction::UserSteppingAction(const G4Step* astep)
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

