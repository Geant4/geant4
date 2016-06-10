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
/// \file runAndEvent/RE01/src/RE01UnknownDecayPhysics.cc
/// \brief Implementation of the RE01UnknownDecayPhysics class
//
// $Id: RE01UnknownDecayPhysics.cc 68761 2013-04-05 12:35:00Z gcosmo $
//
#include "RE01UnknownDecayPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01UnknownDecayPhysics::RE01UnknownDecayPhysics(const G4String& name)
                     :  G4VPhysicsConstructor(name)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01UnknownDecayPhysics::~RE01UnknownDecayPhysics()
{;}

#include "G4ProcessManager.hh"
#include "G4UnknownParticle.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE01UnknownDecayPhysics::ConstructParticle()
{
  G4UnknownParticle::UnknownParticleDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE01UnknownDecayPhysics::ConstructProcess()
{
  // Add Decay Process
  aParticleIterator->reset();
  while( (*aParticleIterator)() ){
    G4ParticleDefinition* particle = aParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if(particle->GetParticleName()=="unknown") {
      pmanager ->AddProcess(&fUnknownDecay);
      pmanager ->SetProcessOrdering(&fUnknownDecay, idxPostStep);
    }
  }
}
