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
// $Id: G4DecayPhysics.cc,v 1.2 2010-06-03 14:37:24 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4DecayPhysics
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard
// 05.12.2005 V.Ivanchenko add controlled verbosity
// 25.04.2006 V.Ivanchenko fix problem of destructor 
//
//----------------------------------------------------------------------------
//

#include "G4DecayPhysics.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"


G4DecayPhysics::G4DecayPhysics(G4int ver)
  :  G4VPhysicsConstructor("Decay"), verbose(ver), wasActivated(false)
{}

G4DecayPhysics::G4DecayPhysics(const G4String& name, G4int ver)
  :  G4VPhysicsConstructor(name), verbose(ver), wasActivated(false)
{}

G4DecayPhysics::~G4DecayPhysics()
{
  if(wasActivated) delete fDecayProcess;
}

void G4DecayPhysics::ConstructParticle()
{

// G4cout << "G4DecayPhysics::ConstructParticle" << G4endl;
  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

void G4DecayPhysics::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  // Add Decay Process
  fDecayProcess = new G4Decay();
  theParticleIterator->reset();
  G4ParticleDefinition* particle=0;
  G4ProcessManager* pmanager=0;

  while( (*theParticleIterator)() )
  {
    particle = theParticleIterator->value();
    pmanager = particle->GetProcessManager();
    if( fDecayProcess->IsApplicable(*particle) ) 
    { 
      if(verbose > 1)
        G4cout << "### Decays for " << particle->GetParticleName() << G4endl;
      pmanager -> AddProcess(fDecayProcess);
      pmanager -> SetProcessOrdering(fDecayProcess, idxPostStep);
      pmanager -> SetProcessOrdering(fDecayProcess, idxAtRest);
    }
  }
}


