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
// $Id: P6DExtDecayerPhysics.cc 100687 2016-10-31 11:20:33Z gcosmo $
//
/// \file eventgenerator/pythia/decayer6/src/P6DExtDecayerPhysics.cc
/// \brief Implementation of the P6DExtDecayerPhysics class
///
/// \author I. Hrivnacova; IPN, Orsay

#include "P6DExtDecayerPhysics.hh"
#include "G4Pythia6Decayer.hh"

#include <G4ParticleDefinition.hh>
#include <G4ProcessManager.hh>
#include <G4Decay.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

P6DExtDecayerPhysics::P6DExtDecayerPhysics(const G4String& name)
  : G4VPhysicsConstructor(name)
{
/// Standard constructor
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

P6DExtDecayerPhysics::~P6DExtDecayerPhysics() 
{
/// Destructor
}

//
// protected methods
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void P6DExtDecayerPhysics::ConstructParticle()
{
/// Nothing to be done here
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void P6DExtDecayerPhysics::ConstructProcess()
{
/// Loop over all particles instantiated and add external decayer
/// to all decay processes if External decayer is set

  // Create Geant4 external decayer
  G4Pythia6Decayer* extDecayer = new G4Pythia6Decayer();
  extDecayer->SetVerboseLevel(1); 
     // The extDecayer will be deleted in G4Decay destructor

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)())
  {    
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    
    if ( verboseLevel > 1 ) {
      G4cout << "Setting ext decayer for: " 
             <<  particleIterator->value()->GetParticleName()
             << G4endl;
    } 
    
    G4ProcessVector* processVector = pmanager->GetProcessList();
    for (G4int i=0; i<processVector->length(); i++) {
    
      G4Decay* decay = dynamic_cast<G4Decay*>((*processVector)[i]);
      if ( decay ) decay->SetExtDecayer(extDecayer);
    }              
  }

  if ( verboseLevel > 0 ) {
    G4cout << "External decayer physics constructed." << G4endl;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
