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
// $Id: G4HadronInelasticQBBC.cc,v 1.18 2009-10-04 16:06:19 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronInelasticQBBC
//
// Author: 2 October 2009 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4HadronInelasticQBBC.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4HadronicInteraction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4UInelasticCrossSection.hh"

#include "G4QGSBuilder.hh"
#include "G4FTFBuilder.hh"

#include "G4QStringChipsParticleLevelInterface.hh"
#include "G4CascadeInterface.hh"
#include "G4BinaryCascade.hh"
#include "G4LCapture.hh"

#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPFission.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronHPFissionData.hh"
#include "G4NeutronHPCaptureData.hh"


G4HadronInelasticQBBC::G4HadronInelasticQBBC(G4int ver, const G4String& type)
  : G4VHadronPhysics("hInelastic"),
    htype(type),verbose(ver),wasActivated(false)
{}

G4HadronInelasticQBBC::G4HadronInelasticQBBC(const G4String&, G4int ver, 
					     G4bool, G4bool,
					     G4bool, G4bool hp, G4bool)
  : G4VHadronPhysics("hInelastic"),verbose(ver),wasActivated(false)
{
  htype = "QBBC";
  if(hp) htype = "QBBC_HP";
}

G4HadronInelasticQBBC::~G4HadronInelasticQBBC()
{}

void G4HadronInelasticQBBC::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  if(verbose > 1) {
    G4cout << "### HadronInelasticQBBC Construct Process with type <"
	   << htype << ">" << G4endl;
  }

  G4bool flagHP = false;
  if(htype == "QBBC_HP") { flagHP = true; }
  //  G4bool flagBinary = false;

  // configure models
  G4HadronicInteraction* theQGSP = 
    BuildModel(new G4QGSBuilder("QGSP",true,false),9.5*GeV,100.*TeV);
  G4HadronicInteraction* theFTFP = 
    BuildModel(new G4FTFBuilder("FTFP"),4.5*GeV,25.*GeV);
  G4HadronicInteraction* theFTFP1 = 
    BuildModel(new G4FTFBuilder("FTFP"),4.5*GeV,100.*TeV);
  G4HadronicInteraction* theBERT = 
    NewModel(new G4CascadeInterface(),0.0,9.9*GeV);
  G4HadronicInteraction* theCHIPS = 
    NewModel(new G4QStringChipsParticleLevelInterface(),0.0,7.5*GeV);

  // loop over particles
  theParticleIterator->reset();
  while( (*theParticleIterator)() ) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String pname = particle->GetParticleName();
    if(verbose > 1) { G4cout << "### HadronInelasticQBBC:  " << pname << G4endl; }
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
       pname == "neutron"   || 
       pname == "omega-"    || 
       pname == "pi-"       || 
       pname == "pi+"       || 
       pname == "proton"    || 
       pname == "sigma-"    || 
       pname == "sigma+"    || 
       pname == "sigma0"    || 
       pname == "xi-"       || 
       pname == "xi0") {
      
      G4HadronicProcess* hp = FindInelasticProcess(particle);
      if(!hp) continue;
      G4ProcessManager* pmanager = particle->GetProcessManager();
      pmanager->AddDiscreteProcess(hp);

      if(pname == "proton") {
	hp->AddDataSet(new G4BGGNucleonInelasticXS(particle));

	hp->RegisterMe(theQGSP);
	//        hp->RegisterMe(theFTFP);
	hp->RegisterMe(theBERT);

      } else if(pname == "neutron") {
	hp->AddDataSet(new G4BGGNucleonInelasticXS(particle));

	hp->RegisterMe(theQGSP);
	//        hp->RegisterMe(theFTFP);
       
	G4HadronicProcess* capture = FindCaptureProcess();
	pmanager->AddDiscreteProcess(capture);

	if(flagHP) {
	  G4HadronicInteraction* theBERT1 = 
	    NewModel(new G4CascadeInterface(), 19.5*MeV, 7.*GeV);
	  hp->RegisterMe(theBERT1);
	  hp->RegisterMe(new G4NeutronHPInelastic());
	  hp->AddDataSet(new G4NeutronHPInelasticData());

	  capture->RegisterMe(new G4NeutronHPCapture());
	  capture->AddDataSet(new G4NeutronHPCaptureData());

	  G4HadronicProcess* fission = FindFissionProcess();
	  pmanager->AddDiscreteProcess(fission);
	  fission->RegisterMe(new G4NeutronHPFission());
	  fission->AddDataSet(new G4NeutronHPFissionData());

	} else {
	  hp->RegisterMe(theBERT);
	  capture->RegisterMe(new G4LCapture());
	}

      } else if(pname == "pi-" || pname == "pi+") {
	hp->AddDataSet(new G4BGGPionInelasticXS(particle));
	hp->RegisterMe(theQGSP);
	// hp->RegisterMe(theFTFP);
        hp->RegisterMe(theBERT);

      } else if(pname == "kaon-"     || 
		pname == "kaon+"     || 
		pname == "kaon0S"    || 
		pname == "kaon0L"    ||
		pname == "lambda"    || 
		pname == "sigma-"    || 
		pname == "sigma+"    || 
		pname == "sigma0"    || 
		pname == "xi-"       || 
		pname == "xi0") {
        hp->RegisterMe(theFTFP1);
        hp->RegisterMe(theBERT);

      } else {
	hp->RegisterMe(theFTFP1);
        hp->RegisterMe(theCHIPS);

      }

      if(verbose > 1) {
	G4cout << "### HadronInelasticQBBC: " << hp->GetProcessName()
	       << " added for " << pname << G4endl;
      }
    }
  }
}
