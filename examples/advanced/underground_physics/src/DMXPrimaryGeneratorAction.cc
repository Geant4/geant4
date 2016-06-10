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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
//
// PrimaryGeneratorAction program
// --------------------------------------------------------------

#include "DMXPrimaryGeneratorAction.hh"

#ifdef DMXENV_GPS_USE
#include "G4GeneralParticleSource.hh"
#else
#include "DMXParticleSource.hh"
#endif

#include "DMXAnalysisManager.hh"

#include "G4Event.hh"

#include "Randomize.hh"

#include "globals.hh"


DMXPrimaryGeneratorAction::DMXPrimaryGeneratorAction() {
  
#ifdef DMXENV_GPS_USE
  particleGun = new G4GeneralParticleSource();
#else
  particleGun = new DMXParticleSource();
#endif

  energy_pri=0;
  //  seeds=NULL;
  seeds[0] =-1;
  seeds[1] =-1;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXPrimaryGeneratorAction::~DMXPrimaryGeneratorAction() {

  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  energy_pri = 0.;

  // seeds
  seeds[0] = *G4Random::getTheSeeds();
  seeds[1] = *(G4Random::getTheSeeds()+1);

  particleGun->GeneratePrimaryVertex(anEvent);

  energy_pri = particleGun->GetParticleEnergy();

  //Fill ntuple #1
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillNtupleDColumn(1,0,energy_pri);
  man->AddNtupleRow(1);
}


