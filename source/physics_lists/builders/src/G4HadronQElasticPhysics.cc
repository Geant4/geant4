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
// $Id: G4HadronQElasticPhysics.cc,v 1.9 2010-06-03 14:28:32 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronQElasticPhysics
//
// Author: 17 Nov 2006 V.Ivanchenko
//
// Modified:
// 03.06.2010 V.Ivanchenko cleanup constructors and ConstructProcess method
//
//----------------------------------------------------------------------------
//
// CHIPS x-sections and generator (G4QElastic) for n and p
// LHEP x-section and generator for the rest 

#include "G4HadronQElasticPhysics.hh"

#include "G4UHadronElasticProcess.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadronElastic.hh"
#include "G4QElastic.hh"

#include "G4VQCrossSection.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

G4HadronQElasticPhysics::G4HadronQElasticPhysics(G4int ver)
  : G4VPhysicsConstructor("hElasticCHIPS_LHEP"), verbose(ver),
    wasActivated(false)
{
  if(verbose > 1) { 
    G4cout << "### G4HadronQElasticPhysics: " << GetPhysicsName() 
	   << G4endl; 
  }
  model = 0;
}

G4HadronQElasticPhysics::G4HadronQElasticPhysics(const G4String&,  G4int ver)
  : G4VPhysicsConstructor("hElasticCHIPS_UELAST"), verbose(ver),
    wasActivated(false)
{
  if(verbose > 1) { 
    G4cout << "### G4HadronQElasticPhysics: " << GetPhysicsName() 
	   << G4endl; 
  }
  model = 0;
}

G4HadronQElasticPhysics::~G4HadronQElasticPhysics()
{
  delete model;
}

void G4HadronQElasticPhysics::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void G4HadronQElasticPhysics::ConstructProcess()
{
  if(wasActivated) { return; }
  wasActivated = true;

  G4double elimit = DBL_MAX;

  if(verbose > 1) {
    G4cout << "### HadronQElasticPhysics: use HE limit " << elimit << " MeV" 
	   << G4endl;
  }
  process = new G4QElastic();

  model = new G4HadronElastic();
  model->SetHEModelLowLimit(elimit);
  G4VQCrossSection* man  = model->GetCS();

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String pname = particle->GetParticleName();
    if(pname == "anti_lambda"  ||
       pname == "anti_neutron" ||
       pname == "anti_omega-"  || 
       pname == "anti_proton"  || 
       pname == "anti_sigma-"  || 
       pname == "anti_sigma+"  || 
       pname == "anti_xi-"  || 
       pname == "anti_xi0"  || 
       pname == "kaon-"     || 
       pname == "kaon+"     || 
       pname == "kaon0S"    || 
       pname == "kaon0L"    || 
       pname == "lambda"    || 
       pname == "omega-"    || 
       pname == "pi-"       || 
       pname == "pi+"       || 
       pname == "sigma-"    || 
       pname == "sigma+"    || 
       pname == "xi-"       || 
       pname == "alpha"     ||
       pname == "deuteron"  ||
       pname == "triton") {
      
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4UHadronElasticProcess* hel = new G4UHadronElasticProcess("hElastic");
      hel->SetQElasticCrossSection(man);
      hel->RegisterMe(model);
      pmanager->AddDiscreteProcess(hel);

    } else if(pname == "neutron" || pname == "proton") {   

      G4ProcessManager* pmanager = particle->GetProcessManager();
      pmanager->AddDiscreteProcess(process);

      if(verbose > 0)
	G4cout << "### QElastic added for " 
	       << particle->GetParticleName() << G4endl;
    }
  }
}


