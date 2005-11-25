//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: GeneralPhysics.cc,v 1.2 2005-11-25 15:38:50 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   GeneralPhysics
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 09.11.2005 G.Folger: don't  keep processes as data members, but new these
//
//----------------------------------------------------------------------------
//
#include "GeneralPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
// 
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"

GeneralPhysics::GeneralPhysics(const G4String& name)
                     :  G4VPhysicsConstructor(name), 
		        wasActivated(false),
			pDecayProcess(0)
{
}

GeneralPhysics::~GeneralPhysics()
{
  if(wasActivated)
  {
    theParticleIterator->reset();
    while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      if (pDecayProcess->IsApplicable(*particle) && pmanager) pmanager ->RemoveProcess(pDecayProcess);
    }
  }
}

void GeneralPhysics::ConstructParticle()
{

 G4cout << "GeneralPhysics::ConstructParticle" << G4endl;
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
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();  
}

void GeneralPhysics::ConstructProcess()
{
  // Add Decay Process
  if (!pDecayProcess ) pDecayProcess=new G4Decay(); 
  theParticleIterator->reset();
  G4ParticleDefinition* particle=0;
  G4ProcessManager* pmanager=0;
  while( (*theParticleIterator)() )
  {
    particle = theParticleIterator->value();
    pmanager = particle->GetProcessManager();
    if( pDecayProcess->IsApplicable(*particle) ) 
    { 
      pmanager -> AddProcess(pDecayProcess);
      pmanager -> SetProcessOrdering(pDecayProcess, idxPostStep);
      pmanager -> SetProcessOrdering(pDecayProcess, idxAtRest);
    }
  }
}


// 2002 by J.P. Wellisch

