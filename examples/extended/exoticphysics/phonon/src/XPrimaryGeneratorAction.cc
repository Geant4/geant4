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
/// \file exoticphysics/phonon/src/XPrimaryGeneratorAction.cc
/// \brief Implementation of the XPrimaryGeneratorAction class
//
// $Id$
//

#include "XPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4RandomDirection.hh"

#include "XTPhononFast.hh"
#include "XTPhononSlow.hh"
#include "XLPhonon.hh"

#include "G4SystemOfUnits.hh"

using namespace std;

XPrimaryGeneratorAction::XPrimaryGeneratorAction()
{ 
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);   

  // default particle kinematic
  particleGun->SetParticleDefinition(XLPhonon::PhononDefinition());
  particleGun->SetParticleMomentumDirection(G4RandomDirection());
  particleGun->SetParticlePosition(G4ThreeVector(0.0,0.0,0.0));
  particleGun->SetParticleEnergy(1e-4*eV);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XPrimaryGeneratorAction::~XPrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 
void XPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
 
  particleGun->SetParticleMomentumDirection(G4RandomDirection());

  G4double selector = G4UniformRand();
  if(selector<0.53539) {
         particleGun->SetParticleDefinition(XTPhononSlow::PhononDefinition()); 
       }
       else if(selector<0.90217) {
         particleGun->SetParticleDefinition(XTPhononFast::PhononDefinition());
       }
       else {
         particleGun->SetParticleDefinition(XLPhonon::PhononDefinition());
  }

  //Set phonon energy.
  //Do not set momentum direction here.
  //Any momentum direction set here will be overwritten
  //by XPhononStackingAction::ClassifyNewTrack
  particleGun->SetParticleEnergy(0.0075*eV);
  particleGun->GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


