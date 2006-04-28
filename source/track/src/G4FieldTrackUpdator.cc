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
// $Id: G4FieldTrackUpdator.cc,v 1.1 2006-04-28 13:45:48 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//   M. Asai - first implementation Apr/28/2006
//
//---------------------------------------------------------------
//
// G4FieldTrackUpdator.cc
//
//---------------------------------------------------------------

#include "globals.hh"
#include "G4FieldTrackUpdator.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4TrackStatus.hh"
#include "G4FieldTrack.hh"

G4FieldTrack* G4FieldTrackUpdator::CreateFieldTrack(const G4Track* trk)
{
  G4FieldTrack* ftrk = new G4FieldTrack(
    trk->GetPosition(),
    trk->GetGlobalTime(),
    trk->GetMomentumDirection(),
    trk->GetKineticEnergy(),
    trk->GetDynamicParticle()->GetMass(),
    trk->GetDynamicParticle()->GetCharge(),
    trk->GetDynamicParticle()->GetPolarization(),
    0.0                   // magnetic dipole moment to be implemented
    );
  return ftrk;
}

void G4FieldTrackUpdator::Update(G4FieldTrack* ftrk,const G4Track* trk)
{
  ftrk->UpdateState(
    trk->GetPosition(),     
    trk->GetGlobalTime(),
    trk->GetMomentumDirection(),
    trk->GetKineticEnergy()
    );
  ftrk->SetChargeAndMoments(
    trk->GetDynamicParticle()->GetCharge()
    );
  ftrk->SetSpin(
    trk->GetDynamicParticle()->GetPolarization()
    );
}


