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
// $Id: G4UnknownDecayPhysics.cc 70999 2013-06-09 01:37:53Z adotti $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4UnknownDecayPhysics
//
// Author: 2016 - M. Asai
//
//----------------------------------------------------------------------------
//

#include "G4UnknownDecayPhysics.hh"

#include "G4ProcessManager.hh"
#include "G4UnknownParticle.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4UnknownDecayPhysics);

G4ThreadLocal G4UnknownDecay* G4UnknownDecayPhysics::fDecayProcess = 0;
G4ThreadLocal G4bool G4UnknownDecayPhysics::wasActivated = false;

G4UnknownDecayPhysics::G4UnknownDecayPhysics(G4int ver)
  :  G4VPhysicsConstructor("UnknownDecay"), verbose(ver)
{
}

G4UnknownDecayPhysics::G4UnknownDecayPhysics(const G4String& name, G4int ver)
  :  G4VPhysicsConstructor(name), verbose(ver)
{
}

G4UnknownDecayPhysics::~G4UnknownDecayPhysics()
{
}

void G4UnknownDecayPhysics::ConstructParticle()
{
  G4UnknownParticle::UnknownParticleDefinition();
}

void G4UnknownDecayPhysics::ConstructProcess()
{
  if(wasActivated) { return; }
  wasActivated = true;

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // Add Decay Process
  fDecayProcess = new G4UnknownDecay();
  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  G4ParticleDefinition* particle=0;

  while( (*myParticleIterator)() )
  {
    particle = myParticleIterator->value();
    if( fDecayProcess->IsApplicable(*particle) ) 
    { 
      if(verbose > 1) {
        G4cout << "### Decays for " << particle->GetParticleName() << G4endl;
      }
      ph->RegisterProcess(fDecayProcess, particle);
    }
  }
}
