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
// $Id: G4HadronElasticPhysics.cc,v 1.11 2010-06-03 14:28:32 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronElasticPhysics
//
// Author: 11 April 2006 V. Ivanchenko
//
// Modified:
// 05.07.2006 V.Ivanchenko define process by particle name; 
//                         fix problem of initialisation of HP
// 24.07.2006 V.Ivanchenko add G4NeutronHPElasticData 
// 10.08.2006 V.Ivanchenko separate neutrons from other particles
// 17.11.2006 V.Ivanchenko do not redefine G4HadronElastic default parameters
// 19.02.2007 V.Ivanchenko set QModelLowLimit and LowestEnergyLimit to zero
// 19.02.2007 A.Howard set QModelLowLimit and LowestEnergyLimit to zero 
//                     for neutrons
// 06.03.2007 V.Ivanchenko use updated interface to G4UElasticCrossSection
// 03.06.2010 V.Ivanchenko cleanup constructors and ConstructProcess method
//
//----------------------------------------------------------------------------
//
// CHIPS for sampling scattering for p and n
// Glauber model for samplimg of high energy pi+- (E > 1GeV)
// LHEP sampling model for the other particle
// CHIPS cross sections for p and n off Hydrogen and Helium targets 
// LHEP cross sections for other particles

#include "G4HadronElasticPhysics.hh"

#include "G4HadronicProcess.hh"
#include "G4HadronicInteraction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4Neutron.hh"

#include "G4UHadronElasticProcess.hh"
#include "G4HadronElastic.hh"

#include "G4VQCrossSection.hh"

G4HadronElasticPhysics::G4HadronElasticPhysics(G4int ver)
  : G4VPhysicsConstructor("hElasticUEL_CHIPS"), verbose(ver),  
    wasActivated(false)
{
  if(verbose > 1) { 
    G4cout << "### G4HadronElasticPhysicsHP: " << GetPhysicsName() 
	   << G4endl; 
  }
}

G4HadronElasticPhysics::G4HadronElasticPhysics(
    const G4String&,  G4int ver, G4bool, G4bool)
  : G4VPhysicsConstructor("hElasticUEL_CHIPS"), verbose(ver),  
    wasActivated(false)
{
  if(verbose > 1) { 
    G4cout << "### G4HadronElasticPhysicsHP: " << GetPhysicsName() 
	   << G4endl; 
  }
}

G4HadronElasticPhysics::~G4HadronElasticPhysics()
{}

void G4HadronElasticPhysics::ConstructParticle()
{
// G4cout << "G4HadronElasticPhysics::ConstructParticle" << G4endl;
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void G4HadronElasticPhysics::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  if(verbose > 1) {
    G4cout << "### HadronElasticPhysics Construct Processes with the model <" 
	   << mname << ">" << G4endl;
  }

  G4HadronElastic* model = new G4HadronElastic();
  G4VQCrossSection* man = model->GetCS();
  model->SetQModelLowLimit(0.0);
  model->SetLowestEnergyLimit(0.0);

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
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
       pname == "sigma-"    || 
       pname == "sigma+"    || 
       pname == "xi-"       || 
       pname == "neutron"   || 
       pname == "proton"    || 
       pname == "pi+"       || 
       pname == "pi-"       || 
       pname == "alpha"     ||
       pname == "deuteron"  ||
       pname == "triton") {
      
      G4UHadronElasticProcess* hel = new G4UHadronElasticProcess("hElastic");
      hel->SetQElasticCrossSection(man);
      hel->RegisterMe(model);
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronElasticPhysics added for " 
	       << particle->GetParticleName() << G4endl;
      }
    }
  }
}


