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
// G4FieldTrackUpdator class implementation
//
// Author: M.Asai, 28 April 2006
// --------------------------------------------------------------------

#include "globals.hh"
#include "G4FieldTrackUpdator.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4TrackStatus.hh"
#include "G4FieldTrack.hh"

//---------------------------------------------------------------------
G4FieldTrack* G4FieldTrackUpdator::CreateFieldTrack(const G4Track* trk)
{
  return new G4FieldTrack(
    trk->GetPosition(), trk->GetGlobalTime(), trk->GetMomentumDirection(),
    trk->GetKineticEnergy(), trk->GetDynamicParticle()->GetMass(),
    trk->GetDynamicParticle()->GetCharge(),
    trk->GetDynamicParticle()->GetPolarization(),
    0.0  // magnetic dipole moment to be implemented
  );
}

//---------------------------------------------------------------------
void G4FieldTrackUpdator::Update(G4FieldTrack* ftrk, const G4Track* trk)
{
  const G4DynamicParticle* ptDynamicParticle = trk->GetDynamicParticle();

  // The following properties must be updated 1) for each new track, and
  ftrk->SetRestMass(ptDynamicParticle->GetMass());
  // 2) Since ion can lose/gain electrons, this must be done at every step

  ftrk->UpdateState(trk->GetPosition(), trk->GetGlobalTime(),
                    trk->GetMomentumDirection(), trk->GetKineticEnergy());

#ifdef G4CHECK
  if((trk->GetMomentum() - ftrk->GetMomentum()).mag2() >
     1.e-16 * trk->GetMomentum().mag2())
  {
    G4cerr << "ERROR> G4FieldTrackUpdator sees *Disagreement* in momentum "
           << G4endl;
    G4cout << "  FTupdator: Tracking Momentum= " << trk->GetMomentum()
           << G4endl;
    G4cout << "  FTupdator: FldTrack Momentum= " << ftrk->GetMomentum()
           << G4endl;
    G4cout << "  FTupdator: FldTrack-Tracking= "
           << ftrk->GetMomentum() - trk->GetMomentum() << G4endl;
  }
#endif

  ftrk->SetProperTimeOfFlight(trk->GetProperTime());

  ftrk->SetChargeAndMoments(ptDynamicParticle->GetCharge(),
                            ptDynamicParticle->GetMagneticMoment());
  ftrk->SetPDGSpin(ptDynamicParticle->GetParticleDefinition()->GetPDGSpin());
  // The charge can change during tracking
  ftrk->SetSpin(ptDynamicParticle->GetPolarization());
}
