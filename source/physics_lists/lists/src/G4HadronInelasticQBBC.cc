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
// $Id: G4HadronInelasticQBBC.cc,v 1.24 2009/11/25 18:55:56 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-03 $
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
//#include "G4UInelasticCrossSection.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4PiNuclearCrossSection.hh"

#include "G4QGSBuilder.hh"
#include "G4FTFBuilder.hh"

#include "G4QStringChipsParticleLevelInterface.hh"
#include "G4CascadeInterface.hh"
#include "G4BinaryCascade.hh"
#include "G4LCapture.hh"
#include "G4NeutronRadCapture.hh"

enum QBBCType
{
  fQBBC = 0,
  fQBBC_XGG,
  fQBBC_XGGSN
};

G4HadronInelasticQBBC::G4HadronInelasticQBBC(const G4String& name, G4int ver, 
    G4bool, G4bool,G4bool, G4bool, G4bool)
  : G4VHadronPhysics("hInelastic"),verbose(ver),wasActivated(false)
{
  htype = name;
}

G4HadronInelasticQBBC::~G4HadronInelasticQBBC()
{}

void G4HadronInelasticQBBC::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  // select type
  QBBCType fType = fQBBC;
  if("QBBC_XGG" == htype)        { fType = fQBBC_XGG; }
  else if("QBBC_XGGSN" == htype) { fType = fQBBC_XGGSN; }

  if(verbose > 1) {
    G4cout << "### HadronInelasticQBBC Construct Process with type <"
	   << htype << ">" << G4endl;
  }

  // configure models
  G4HadronicInteraction* theQGSP = 
    BuildModel(new G4QGSBuilder("QGSP",true,false),9.5*GeV,100.*TeV);
  G4HadronicInteraction* theFTFP = 
    BuildModel(new G4FTFBuilder("FTFP"),4.5*GeV,25.*GeV);
  G4HadronicInteraction* theFTFP1 = 
    BuildModel(new G4FTFBuilder("FTFP"),4.5*GeV,100.*TeV);
  G4HadronicInteraction* theBERT = 
    NewModel(new G4CascadeInterface(),0.0,6.5*GeV);
  //G4HadronicInteraction* theBERT1 = 
  //  NewModel(new G4CascadeInterface(),2.5*GeV,6.5*GeV);
  //G4HadronicInteraction* theBIC = 
  //  NewModel(new G4BinaryCascade(),0.0,3.5*GeV);
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
      
      // common for all particles
      G4HadronicProcess* hp = FindInelasticProcess(particle);
      if(!hp) { continue; }
      G4ProcessManager* pmanager = particle->GetProcessManager();
      pmanager->AddDiscreteProcess(hp);

      // model and X-section configuration
      if(pname == "proton") {
        if(fType == fQBBC) { 
	  hp->AddDataSet(new G4ProtonInelasticCrossSection()); 
	} else {
          hp->AddDataSet(new G4BGGNucleonInelasticXS(particle));
	}
	hp->RegisterMe(theQGSP);
	hp->RegisterMe(theFTFP);
	hp->RegisterMe(theBERT);
	//	hp->RegisterMe(theBERT1);
	//	hp->RegisterMe(theBIC);

      } else if(pname == "neutron") {
        if(fType == fQBBC) { 
	  hp->AddDataSet(new G4NeutronInelasticCrossSection()); 
	} else if(fType == fQBBC_XGG) {
          hp->AddDataSet(new G4BGGNucleonInelasticXS(particle));
	} else {
          hp->AddDataSet(new G4NeutronInelasticXS());
	}

	hp->RegisterMe(theQGSP);
	hp->RegisterMe(theFTFP);
	//	hp->RegisterMe(theBERT1);
       
	G4HadronicProcess* capture = FindCaptureProcess();
	pmanager->AddDiscreteProcess(capture);

	hp->RegisterMe(theBERT);
	//hp->RegisterMe(theBIC);
	capture->RegisterMe(new G4NeutronRadCapture());
	//capture->RegisterMe(new G4LCapture());
	if(fType == fQBBC_XGGSN) {
	  capture->AddDataSet(new G4NeutronCaptureXS());
	}

      } else if(pname == "pi-" || pname == "pi+") {
        if(fType == fQBBC) { 
	  hp->AddDataSet(new G4PiNuclearCrossSection()); 
	} else {
          hp->AddDataSet(new G4BGGPionInelasticXS(particle));
	}
	hp->RegisterMe(theQGSP);
	hp->RegisterMe(theFTFP);
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
