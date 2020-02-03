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
/// \file eventgenerator/HepMC/HepMCEx02/src/H02PrimaryGeneratorAction.cc
/// \brief Implementation of the H02PrimaryGeneratorAction class
//
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "HepMCG4AsciiReader.hh"
#include "HepMCG4PythiaInterface.hh"
#include "H02PrimaryGeneratorAction.hh"
#include "H02PrimaryGeneratorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02PrimaryGeneratorAction::H02PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction()
{
  // default generator is particle gun.
  fCurrentGenerator= fParticleGun= new G4ParticleGun();
  fCurrentGeneratorName= "particleGun";
  fHepmcAscii= new HepMCG4AsciiReader();
#ifdef G4LIB_USE_PYTHIA
  fPythiaGen= new HepMCG4PythiaInterface();
#else
  fPythiaGen= 0;
#endif

  fGentypeMap["particleGun"]= fParticleGun;
  fGentypeMap["hepmcAscii"]= fHepmcAscii;
  fGentypeMap["pythia"]= fPythiaGen;

  fMessenger= new H02PrimaryGeneratorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02PrimaryGeneratorAction::~H02PrimaryGeneratorAction()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H02PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(fCurrentGenerator)
    fCurrentGenerator-> GeneratePrimaryVertex(anEvent);
  else
    G4Exception("H02PrimaryGeneratorAction::GeneratePrimaries",
                "InvalidSetup", FatalException,
                "Generator is not instanciated.");
}
