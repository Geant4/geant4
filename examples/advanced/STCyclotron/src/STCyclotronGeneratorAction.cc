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
//  Author: F. Poignant, floriane.poignant@gmail.com
// 
// file STyclotronPrimaryGeneratorAction.cc

#include "STCyclotronPrimaryGeneratorAction.hh"
#include "STCyclotronRun.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "STCyclotronPrimaryGeneratorActionMessenger.hh"

STCyclotronPrimaryGeneratorAction::STCyclotronPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction()
{
  fMessenger = new STCyclotronPrimaryGeneratorActionMessenger(this);
  fParticleBeam  = new G4GeneralParticleSource();
  fBeamCurrent = 10.E-6 ; //ampere;
  //The rest of the parameters (type of particle, energy, beam shape ..) are defined in the init_beam.vis class.
}

STCyclotronPrimaryGeneratorAction::~STCyclotronPrimaryGeneratorAction()
{
  delete fMessenger;
  delete fParticleBeam;
}

void STCyclotronPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //Set up the number of particles per event 
  G4double timePerEvent = 1.E-11 ; //s;
  G4double chargeParticle = fParticleBeam->GetParticleDefinition()->GetPDGCharge()*1.6E-19;
  G4double numberOfPart = std::abs(fBeamCurrent*timePerEvent/chargeParticle); 
  G4String name = fParticleBeam->GetParticleDefinition()->GetParticleName();
  G4double energy = fParticleBeam->GetParticleEnergy();

  G4int fPrimariesPerEvent = (G4int)numberOfPart;
  
  if(fPrimariesPerEvent < 1){
    G4cout << "Warning: number of particles per event below 0: " << numberOfPart << G4endl;
    return;
  }

  fParticleBeam->SetNumberOfParticles(fPrimariesPerEvent);
  fParticleBeam->GeneratePrimaryVertex(anEvent);
  
  STCyclotronRun* fRun = static_cast<STCyclotronRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  fRun->SetPrimariesPerEvent(fPrimariesPerEvent);
  fRun->SetTimePerEvent(timePerEvent);

  fRun->SetBeamName(name);
  fRun->SetBeamCurrent(fBeamCurrent);
  fRun->SetBeamEnergy(energy);
  
  //G4cout << "The new beam current is the following : " << fBeamCurrent << " Ampere." << G4endl;
  //G4cout << "Particles per event : " << numberOfParticlePerEvent << " particles." << G4endl;

}

void STCyclotronPrimaryGeneratorAction::SetBeamCurrent(G4double current)
{
  if(fBeamCurrent!=current){
    fBeamCurrent=current;
    G4cout << "The new beam current is the following : " << fBeamCurrent << " Ampere." << G4endl;
  } 
}
