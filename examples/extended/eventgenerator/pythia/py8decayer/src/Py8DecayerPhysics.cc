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
/// \file eventgenerator/pythia/pythia8decayer/src/Py8DecayerPhysics.cc
/// \brief Implementation of the Py8DecayerPhysics class
///
/// \author J. Yarba; FNAL

#include "Py8DecayerPhysics.hh"
#include "Py8Decayer.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Decay.hh"
#include "G4DecayTable.hh"

// factory
//
#include "G4PhysicsConstructorFactory.hh"
//
// register it with contructor factory
//
G4_DECLARE_PHYSCONSTR_FACTORY(Py8DecayerPhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Py8DecayerPhysics::Py8DecayerPhysics(G4int)
  : G4VPhysicsConstructor("Py8DecayerPhysics")
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Py8DecayerPhysics::~Py8DecayerPhysics() 
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Py8DecayerPhysics::ConstructParticle()
{
   // Nothing needs to be done here
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Py8DecayerPhysics::ConstructProcess()
{
   // Adding external decayer to G4Decay process (per each thread).
   // G4Decay will use the external decayer if G4Decay process is
   // assigned to an unstable particle and that particle does not
   // have its decay table.

   // Loop over all particles instantiated and remove already-assigned
   // decay table for tau's and B+/- so that they will decay through
   // the external decayer (Pythia8).

   // NOTE: The extDecayer will be deleted in G4Decay destructor
   
   Py8Decayer* extDecayer = new Py8Decayer();
   G4bool setOnce = true;

   auto particleIterator=GetParticleIterator();
   particleIterator->reset();
   while ((*particleIterator)())
   {    
      G4ParticleDefinition* particle = particleIterator->value();

      // remove native/existing decay table for
      // a)tau's 
      // b) B+/- 
      // so that G4Decay will use the external decayer
      if ( std::abs(particle->GetPDGEncoding()) == 15 ||
      std::abs(particle->GetPDGEncoding()) == 521 )
      {
        if ( particle->GetDecayTable() )
        {
          delete particle->GetDecayTable();
          particle->SetDecayTable(nullptr);
/*
          if ( verboseLevel > 1 ) {
             G4cout << "Use ext decayer for: " 
                <<  particleIterator->value()->GetParticleName()
                << G4endl;
          } 
*/    
        }
      }

      if(setOnce)
      // One G4Decay object is shared by all unstable particles (per thread).
      // Thus, we set the external decayer only once.
      {
        G4ProcessManager* pmanager = particle->GetProcessManager();    
        G4ProcessVector* processVector = pmanager->GetProcessList();
        for ( size_t i=0; i<processVector->length(); ++i ) 
        {    
           G4Decay* decay = dynamic_cast<G4Decay*>((*processVector)[i]);
           if ( decay ) 
           {
             decay->SetExtDecayer(extDecayer);
             setOnce = false;
           }
        }
      }              
   }

   return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
