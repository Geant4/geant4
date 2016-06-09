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
// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// * LISAPrimaryGeneratorAction class                                 *
// *                                                                  *
// ********************************************************************
//
// HISTORY
// 22/02/2004: migrated from LISA-V04
//
// ********************************************************************


#include "LISAPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"

#include "G4Event.hh"
#include "Randomize.hh"


LISAPrimaryGeneratorAction::LISAPrimaryGeneratorAction() {
  
  particleGun = new G4GeneralParticleSource();

  energy_pri=0;
  seeds[0]=-1;
  seeds[1]=-1;

}


LISAPrimaryGeneratorAction::~LISAPrimaryGeneratorAction() {

  delete particleGun;
}


void LISAPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  energy_pri = 0.;

  // seeds
  seeds[0] = *CLHEP::HepRandom::getTheSeeds();
  seeds[1] = *(CLHEP::HepRandom::getTheSeeds()+1);
  //  G4cout << " 1st seed: " << *seeds << G4endl;;
  //  G4cout << " 2nd seed: " << *(seeds+1) << G4endl;
  // HepRandom::showEngineStatus();
      
  particleGun->GeneratePrimaryVertex(anEvent);

  energy_pri = particleGun->GetParticleEnergy();

}

