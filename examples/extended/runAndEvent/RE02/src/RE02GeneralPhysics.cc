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
// $Id: RE02GeneralPhysics.cc,v 1.2 2006-06-29 17:45:16 gunter Exp $
// --------------------------------------------------------------
//
// 22-Nov-2004 Construt ALL Particles by T. Koi


#include "RE02GeneralPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>

RE02GeneralPhysics::RE02GeneralPhysics(const G4String& name)
                     :  G4VPhysicsConstructor(name)
{
}

RE02GeneralPhysics::~RE02GeneralPhysics()
{
}

#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

void RE02GeneralPhysics::ConstructParticle()
{
   // In Alphabetical Order 
   
   //  Construct all barions
   G4BaryonConstructor* baryonConstructor = new G4BaryonConstructor();
   baryonConstructor -> ConstructParticle();
   delete baryonConstructor;

   // Construct all bosons (including geantinos)
   G4BosonConstructor* bosonConstructor = new G4BosonConstructor();
   bosonConstructor -> ConstructParticle();
   delete bosonConstructor;

   // Construct all ions 
   G4IonConstructor* ionConstructor = new G4IonConstructor();
   ionConstructor -> ConstructParticle();
   delete ionConstructor;

   // Construct all leptons 
   G4LeptonConstructor* leptonConstructor = new G4LeptonConstructor();
   leptonConstructor -> ConstructParticle();
   delete leptonConstructor;

   // Construct all mesons
   G4MesonConstructor* mesonConstructor = new G4MesonConstructor();
   mesonConstructor -> ConstructParticle();
   delete mesonConstructor;

   //  Construct  resonaces and quarks
   G4ShortLivedConstructor* shortLivedConstructor = new G4ShortLivedConstructor();
   shortLivedConstructor -> ConstructParticle();
   delete shortLivedConstructor;

}

#include "G4Decay.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

void RE02GeneralPhysics::ConstructProcess()
{
  // Add Decay Process
   G4Decay* theDecayProcess = new G4Decay();  
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) {
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}


