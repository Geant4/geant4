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
/// \file exoticphysics/phonon/src/XPhononStackingAction.cc
/// \brief Implementation of the XPhononStackingAction class
//
// $Id$
//

#include "XPhononStackingAction.hh"
#include "XPhononTrackInformation.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "XLatticeManager3.hh"
#include "G4RandomDirection.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "XTPhononFast.hh"
#include "XTPhononSlow.hh"
#include "XLPhonon.hh"

XPhononStackingAction::XPhononStackingAction()
{;}

XPhononStackingAction::~XPhononStackingAction()
{;}

G4ClassificationOfNewTrack XPhononStackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
  //This stacking action is necessary to ensure that velocity and 
  //propagation direction are set properly for phonons created with
  //G4ParticleGun

  G4ClassificationOfNewTrack classification = fUrgent;

  if(aTrack->GetParentID() == 0)
    {
      //Obtain LatticeManager for phonon dynamics
      XLatticeManager3* LM = XLatticeManager3::GetXLatticeManager();
      
      //cast to non-const pointer so we can set the velocity
      G4Track* theTrack=(G4Track*)aTrack;

      G4int pol=0;
      if(theTrack->GetDefinition()==XLPhonon::PhononDefinition()) {
        pol = 0;
      } else if(theTrack->GetDefinition()==XTPhononSlow::PhononDefinition()) {
        pol = 1;
      } else if(theTrack->GetDefinition()==XTPhononFast::PhononDefinition()) {
        pol = 2;
      }
      
      //Compute random wave-vector
      G4ThreeVector Ran = G4RandomDirection();

      //Store wave-vector as track information
      theTrack->SetUserInformation(new XPhononTrackInformation(Ran));     
      
      //Compute direction of propagation from wave vector
      G4ThreeVector momentumDir;
      momentumDir =  LM->MapKtoVDir(theTrack->GetVolume(),pol,Ran);
      theTrack->SetMomentumDirection(momentumDir);
      
      //Compute true velocity of propagation
      G4double velocity;
      velocity = LM->MapKtoV(theTrack->GetVolume(), pol, Ran)*m/s;
      theTrack->SetVelocity(velocity);
      theTrack->UseGivenVelocity(true);
      
    }
  
  return classification; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

