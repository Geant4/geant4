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
// $Id: G4HadronInelasticQBBC.cc,v 1.29 2010-11-19 19:50:15 vnivanch Exp $
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

#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4QInelastic.hh"
#include "G4HadronicProcessStore.hh"

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

#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4Evaporation.hh"

enum QBBCType
{
  fQBBC = 0,     // default QBBC
  fQBBC_XGG,     // neutron x-sections from BGG
  fQBBC_XGGSN    // n, p, pi+- x-sections from QGSP_BERT
};

G4HadronInelasticQBBC::G4HadronInelasticQBBC(G4int ver) 
  : G4VHadronPhysics("hInelastic"),verbose(ver),wasActivated(false)
{
  htype = "QBBC";
  theHandler = 0;
  theEvaporation = 0;
}

G4HadronInelasticQBBC::G4HadronInelasticQBBC(const G4String& name, G4int ver, 
    G4bool, G4bool,G4bool, G4bool, G4bool)
  : G4VHadronPhysics("hInelastic"),verbose(ver),wasActivated(false)
{
  htype = name;
  theHandler = 0;
  theEvaporation = 0;
}

G4HadronInelasticQBBC::~G4HadronInelasticQBBC()
{
  delete theHandler;
  delete theEvaporation;
}

void G4HadronInelasticQBBC::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  if(verbose > 1) {
    G4cout << "### HadronInelasticQBBC Construct Process with type <"
	   << htype << ">" << G4endl;
  }

  // PreCompound and Evaporation models are instantiated here
  //theEvaporation = new G4Evaporation();
  //theEvaporation->SetCombinedChannel();
  theHandler = new G4ExcitationHandler();
  //theHandler->SetEvaporation(theEvaporation);
  //theHandler->SetMinEForMultiFrag(3.0*GeV);
  //theHandler->SetMaxAandZForFermiBreakUp(17,9);
  G4PreCompoundModel* thePreCompound = new G4PreCompoundModel(theHandler);

  // configure models
  G4HadronicInteraction* theQGSP = 
    BuildModel(new G4QGSBuilder("QGSP",thePreCompound,true,false),12.5*GeV,100.*TeV);
  G4HadronicInteraction* theFTFP = 
    BuildModel(new G4FTFBuilder("FTFP",thePreCompound),4.0*GeV,25.*GeV);
  G4HadronicInteraction* theFTFP1 = 
    BuildModel(new G4FTFBuilder("FTFP",thePreCompound),4.0*GeV,100.*TeV);

  G4HadronicInteraction* theBERT = 
    NewModel(new G4CascadeInterface(),1.0*GeV,5.0*GeV);
  G4HadronicInteraction* theBERT1 = 
    NewModel(new G4CascadeInterface(),0.0*GeV,5.0*GeV);

  G4BinaryCascade* bic = new G4BinaryCascade();
  bic->SetDeExcitation(thePreCompound);
  G4HadronicInteraction* theBIC = NewModel(bic,0.0,1.5*GeV);

  G4QInelastic* theCHIPS = new G4QInelastic();
  G4HadronicProcessStore* store = G4HadronicProcessStore::Instance();
  store->RegisterExtraProcess(theCHIPS);

  // loop over particles
  theParticleIterator->reset();
  while( (*theParticleIterator)() ) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String pname = particle->GetParticleName();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if(verbose > 1) { 
      G4cout << "### HadronInelasticQBBC:  " << pname << G4endl; 
    }

    //
    // model and X-section configuration per particle type
    //
    if(pname == "proton") {
      G4HadronicProcess* hp = FindInelasticProcess(particle);
      hp->AddDataSet(new G4BGGNucleonInelasticXS(particle));
      
      hp->RegisterMe(theQGSP);
      hp->RegisterMe(theFTFP);
      hp->RegisterMe(theBERT);
      hp->RegisterMe(theBIC);

    } else if(pname == "neutron") {
      G4HadronicProcess* hp = FindInelasticProcess(particle);
      hp->AddDataSet(new G4NeutronInelasticXS());
      hp->RegisterMe(theQGSP);
      hp->RegisterMe(theFTFP);
       
      G4HadronicProcess* capture = FindCaptureProcess();
      capture->AddDataSet(new G4NeutronCaptureXS());
      hp->RegisterMe(theBERT);
      hp->RegisterMe(theBIC);
      capture->RegisterMe(new G4NeutronRadCapture());

    } else if(pname == "pi-" || pname == "pi+") {
      G4HadronicProcess* hp = FindInelasticProcess(particle);
      hp->AddDataSet(new G4BGGPionInelasticXS(particle));
      hp->RegisterMe(theQGSP);
      hp->RegisterMe(theFTFP);
      hp->RegisterMe(theBERT1);

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
      G4HadronicProcess* hp = FindInelasticProcess(particle);
      hp->RegisterMe(theFTFP1);
      hp->RegisterMe(theBERT1);

    } else if(pname == "anti_lambda"  ||
              pname == "anti_neutron" ||
              pname == "anti_omega-"  || 
              pname == "anti_proton"  || 
              pname == "anti_sigma-"  || 
              pname == "anti_sigma+"  || 
              pname == "anti_xi-"  || 
              pname == "anti_xi0"  ||
              pname == "omega-") {
      //G4HadronicProcess* hp = FindInelasticProcess(particle);
      //hp->RegisterMe(theFTFP1);
      //hp->RegisterMe(theCHIPS);
      pmanager->AddDiscreteProcess(theCHIPS);
      store->RegisterParticleForExtraProcess(theCHIPS,particle);
    } 
  }
}
