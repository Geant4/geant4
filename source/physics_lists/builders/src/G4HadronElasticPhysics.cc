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
// $Id: G4HadronElasticPhysics.cc,v 1.16 2010-09-24 12:54:13 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronElasticPhysics 
//
// Author: 23 November 2006 V. Ivanchenko
//
// Modified:
// 21.03.2007 V.Ivanchenko Use G4BGGNucleonElasticXS and G4BGGPionElasticXS; 
//                         Reduce thresholds for HE and Q-models to zero
// 03.06.2010 V.Ivanchenko cleanup constructors and ConstructProcess method
// 29.07.2010 V.Ivanchenko rename this class from G4HadronHElasticPhysics to
//                         G4HadronElasticPhysics, old version of the class
//                         is renamed to G4HadronElasticPhysics93
//
//----------------------------------------------------------------------------
//
// CHIPS for sampling scattering for p and n
// Glauber model for samplimg of high energy pi+- (E > 1GeV)
// LHEP sampling model for the other particle
// BBG cross sections for p, n and pi+- 
// LHEP cross sections for other particles

#include "G4HadronElasticPhysics.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4WHadronElasticProcess.hh"
#include "G4HadronElastic.hh"
#include "G4CHIPSElastic.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "G4AntiNuclElastic.hh"

#include "G4BGGNucleonElasticXS.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4NeutronElasticXS.hh"
#include "G4CHIPSElasticXS.hh"

#include "G4ComponentAntiNuclNuclearXS.hh"  
#include "G4CrossSectionElastic.hh"

G4HadronElasticPhysics::G4HadronElasticPhysics(G4int ver)
  : G4VPhysicsConstructor("hElasticWEL_CHIPS"), verbose(ver), 
    wasActivated(false)
{
  if(verbose > 1) { 
    G4cout << "### G4HadronElasticPhysics: " << GetPhysicsName() 
	   << G4endl; 
  }
}

G4HadronElasticPhysics::G4HadronElasticPhysics(const G4String&,
    G4int ver, G4bool, const G4String&)
  : G4VPhysicsConstructor("hElasticWEL_CHIPS"), verbose(ver), 
    wasActivated(false)
{
  if(verbose > 1) { 
    G4cout << "### G4HadronElasticPhysics: " << GetPhysicsName() 
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

  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void G4HadronElasticPhysics::ConstructProcess()
{
  if(wasActivated) { return; }
  wasActivated = true;

  G4double elimitPi = 1.0*GeV;
  G4double elimitAntiNuc = 100*MeV;
  if(verbose > 1) {
    G4cout << "### HadronElasticPhysics Construct Processes with the limit for pi " 
	   << elimitPi/GeV << " GeV" 
	   << "                                                  for anti-neuclei " 
	   << elimitAntiNuc/GeV << " GeV"	   << G4endl;
  }

  G4AntiNuclElastic* anuc = new G4AntiNuclElastic();
  anuc->SetMinEnergy(elimitAntiNuc);
  G4CrossSectionElastic* anucxs = 
    new G4CrossSectionElastic(anuc->GetComponentCrossSection());

  G4HadronElastic* lhep0 = new G4HadronElastic();
  G4HadronElastic* lhep1 = new G4HadronElastic();
  G4HadronElastic* lhep2 = new G4HadronElastic();
  lhep1->SetMaxEnergy(elimitPi);
  lhep2->SetMaxEnergy(elimitAntiNuc);

  G4CHIPSElastic* chipsp = new G4CHIPSElastic();
  neutronModel = new G4CHIPSElastic();

  G4ElasticHadrNucleusHE* he = new G4ElasticHadrNucleusHE(); 
  he->SetMinEnergy(elimitPi);

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String pname = particle->GetParticleName();
    if(pname == "anti_lambda"  ||
       pname == "anti_neutron" ||
       pname == "anti_omega-"  || 
       pname == "anti_sigma-"  || 
       pname == "anti_sigma+"  || 
       pname == "anti_xi-"  || 
       pname == "anti_xi0"  || 
       pname == "lambda"    || 
       pname == "omega-"    || 
       pname == "sigma-"    || 
       pname == "sigma+"    || 
       pname == "xi-"       || 
       pname == "alpha"     ||
       pname == "deuteron"  ||
       pname == "triton"   
       ) {
      
      G4WHadronElasticProcess* hel = new G4WHadronElasticProcess();
      hel->RegisterMe(lhep0);
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if(pname == "proton") {   

      G4WHadronElasticProcess* hel = new G4WHadronElasticProcess();
      //hel->AddDataSet(new G4BGGNucleonElasticXS(particle));
      hel->AddDataSet(new G4CHIPSElasticXS());
      hel->RegisterMe(chipsp);
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if(pname == "neutron") {   

      neutronProcess = new G4WHadronElasticProcess();
      //neutronProcess->AddDataSet(new G4BGGNucleonElasticXS(particle));
      neutronProcess->AddDataSet(new G4CHIPSElasticXS());
      neutronProcess->RegisterMe(neutronModel);
      pmanager->AddDiscreteProcess(neutronProcess);
      if(verbose > 1) {
	G4cout << "### HadronElasticPhysics: " 
	       << neutronProcess->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if (pname == "pi+" || pname == "pi-") { 

      G4WHadronElasticProcess* hel = new G4WHadronElasticProcess();
      hel->AddDataSet(new G4BGGPionElasticXS(particle));
      hel->RegisterMe(lhep1);
      hel->RegisterMe(he);
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if(pname == "kaon-"     || 
	      pname == "kaon+"     || 
	      pname == "kaon0S"    || 
	      pname == "kaon0L" 
	      ) {
      
      G4WHadronElasticProcess* hel = new G4WHadronElasticProcess();
      hel->RegisterMe(lhep0);
      pmanager->AddDiscreteProcess(hel);
      if(verbose > 1) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if(
       pname == "anti_proton"    || 
       pname == "anti_alpha"     ||
       pname == "anti_deuteron"  ||
       pname == "anti_triton"    ||
       pname == "anti_He3"       ) {

      G4WHadronElasticProcess* hel = new G4WHadronElasticProcess();
      hel->AddDataSet(anucxs);
      hel->RegisterMe(lhep2);
      hel->RegisterMe(anuc);
      pmanager->AddDiscreteProcess(hel);
    }
  }
}


