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
// $Id: G4HadronInelasticQBBC.cc,v 1.2 2007/04/16 11:57:40 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronInelasticQBBC
//
// Author: 11 April 2006 V. Ivanchenko
//
// Modified:
// 05.07.2006 V.Ivanchenko fix problem of initialisation of HP
// 15.04.2007 V.Ivanchenko include quasi-elastic and change FTF low energy
//
//----------------------------------------------------------------------------
//

#include "G4HadronInelasticQBBC.hh"

#include "G4HadronInelasticProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"

#include "G4PiNuclearCrossSection.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionInelasticXS.hh"

#include "G4TheoFSGenerator.hh"
#include "G4StringChipsParticleLevelInterface.hh"
#include "G4StringChipsInterface.hh"
#include "G4QGSMFragmentation.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"

#include "G4BinaryCascade.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"

#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPFission.hh"
#include "G4NeutronHPCapture.hh"

#include "G4HadronProcessStore.hh"
#include "G4UInelasticCrossSection.hh"

G4HadronInelasticQBBC::G4HadronInelasticQBBC(const G4String& name, 
    G4int ver, G4bool ftf, G4bool bert, G4bool chips, G4bool hp, G4bool glauber)
  : G4VPhysicsConstructor(name), verbose(ver), ftfFlag(ftf), bertFlag(bert), 
    chipsFlag(chips), hpFlag(hp), glFlag(glauber), wasActivated(false)
{
  if(verbose > -1) G4cout << "### HadronInelasticQBBC" << G4endl;
  store = G4HadronProcessStore::Instance();
  theHPXSecI = 0;
  theHPXSecC = 0;
  theHPXSecF = 0;
  theCascade = 0;
  theCHIPSCascade = 0;
  theQuasiElastic = 0;
  thePreEquilib = 0;
}

G4HadronInelasticQBBC::~G4HadronInelasticQBBC()
{
  if(wasActivated) {
    if(theCascade) delete theCascade;
    if(theCHIPSCascade) delete theCHIPSCascade;
    if(thePreEquilib) delete thePreEquilib;
    if(theQuasiElastic) delete theQuasiElastic;
    delete theQGStringDecay;
    delete theQGStringModel;
    delete theFTFStringDecay;
    delete theFTFStringModel;
    delete theHPXSecI;
    delete theHPXSecC;
    delete theHPXSecF;
  }
}

void G4HadronInelasticQBBC::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();
}

void G4HadronInelasticQBBC::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  if(verbose > 1) 
    G4cout << "### HadronInelasticQBBC Construct Process" << G4endl;

  G4double minEstring  = 11.5*GeV;
  G4double maxEcascade =  4.*GeV;
  G4double minFtf      =  3.5*GeV;
  G4double maxFtf      = 12.*GeV;

  //Bertini
  G4HadronicInteraction* theBERT = new G4CascadeInterface();
  theBERT->SetMinEnergy(0.0);
  theBERT->SetMaxEnergy(maxEcascade);

  //CHIPS
  G4HadronicInteraction* theCHIPS = new G4StringChipsInterface();
  theCHIPS->SetMinEnergy(0.0);
  theCHIPS->SetMaxEnergy(maxEcascade);

  //QGS
  G4TheoFSGenerator* theQGSModel = new G4TheoFSGenerator();
  theQGStringModel  = new G4QGSModel< G4QGSParticipants >;
  theQGStringDecay  = new G4ExcitedStringDecay(new G4QGSMFragmentation());
  theQGStringModel->SetFragmentationModel(theQGStringDecay);

  //QGSP
  theCascade = new G4GeneratorPrecompoundInterface;
  thePreEquilib = new G4PreCompoundModel(new G4ExcitationHandler);
  theCascade->SetDeExcitation(thePreEquilib);  
  theQGSModel->SetTransport(theCascade);

  //QGSC
  //  theCHIPSCascade   = new G4StringChipsParticleLevelInterface;
  //  theQGSModel->SetTransport(theCHIPSCascade);

  //QGS
  theQuasiElastic = new G4QuasiElasticChannel;
  theQGSModel->SetQuasiElasticChannel(theQuasiElastic);
  theQGSModel->SetHighEnergyGenerator(theQGStringModel);
  theQGSModel->SetMinEnergy(minEstring);
  theQGSModel->SetMaxEnergy(100*TeV);

  //FTF
  G4TheoFSGenerator* theFTFModel = new G4TheoFSGenerator();
  theFTFStringModel = new G4FTFModel();
  theFTFStringDecay = new G4ExcitedStringDecay(new G4LundStringFragmentation());
  theFTFStringModel->SetFragmentationModel(theFTFStringDecay);
  //  theFTFModel->SetTransport(theCHIPSCascade);
  theFTFModel->SetTransport(theCascade);
  theFTFModel->SetHighEnergyGenerator(theFTFStringModel);
  theFTFModel->SetMinEnergy(minFtf);
  if(ftfFlag) theFTFModel->SetMaxEnergy(100*TeV);
  else        theFTFModel->SetMaxEnergy(maxFtf);

  G4TheoFSGenerator* theFTFModel1 = new G4TheoFSGenerator();
  theFTFModel1->SetTransport(theCascade);
  theFTFModel1->SetHighEnergyGenerator(theFTFStringModel);
  theFTFModel1->SetMinEnergy(minFtf);
  theFTFModel1->SetMaxEnergy(100*TeV);

  theFTFModel->SetQuasiElasticChannel(theQuasiElastic);
  theFTFModel1->SetQuasiElasticChannel(theQuasiElastic);


  theParticleIterator->reset();
  while( (*theParticleIterator)() ) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String pname = particle->GetParticleName();
    if(verbose > 1) G4cout << "### HadronInelasticQBBC:  " << pname << G4endl;
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
       pname == "xi-"       || 
       pname == "xi0") {
      
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4HadronInelasticProcess* hp = 
	new G4HadronInelasticProcess("hInelastic", particle);
      pmanager->AddDiscreteProcess(hp);

      if(pname == "proton") {
	hp->AddDataSet(&theXSecP);
        if(ftfFlag) Register(particle,hp,theFTFModel,"FTF");
	else {
	  Register(particle,hp,theQGSModel,"QGS");
	  Register(particle,hp,theFTFModel,"FTF");
	}

	if(bertFlag)       Register(particle,hp,theBERT,"Bertini");
        else if(chipsFlag) Register(particle,hp,theCHIPS,"CHIPS");
	else {
	  G4HadronicInteraction* theBIC = new G4BinaryCascade();
	  theBIC->SetMinEnergy(0.0);
	  theBIC->SetMaxEnergy(maxEcascade);
	  Register(particle,hp,theBIC,"Binary");
	}
	if(glFlag) 
	  hp->AddDataSet(new G4BGGNucleonInelasticXS(particle));

      } else if(pname == "neutron") {
	hp->AddDataSet(&theXSecN);
        if(ftfFlag) Register(particle,hp,theFTFModel,"FTF");
	else {
	  Register(particle,hp,theQGSModel,"QGS");
	  Register(particle,hp,theFTFModel,"FTF");
	}

	G4HadronCaptureProcess* theNeutronCapture = 
	  new G4HadronCaptureProcess("nCapture");
	G4HadronFissionProcess* theNeutronFission = 
	  new G4HadronFissionProcess("nFission");
	pmanager->AddDiscreteProcess(theNeutronCapture);
	pmanager->AddDiscreteProcess(theNeutronFission);

	G4double emin = 0.0;
	if(hpFlag) {
	  emin = 19.5*MeV;
          theHPXSecI = new G4NeutronHPInelasticData;
          theHPXSecC = new G4NeutronHPCaptureData;
	  theHPXSecF = new G4NeutronHPFissionData;
	  hp->AddDataSet(theHPXSecI);
	  theNeutronCapture->AddDataSet(theHPXSecC);
	  theNeutronFission->AddDataSet(theHPXSecF);
          G4NeutronHPInelastic* hpi = new G4NeutronHPInelastic();
          G4NeutronHPCapture* hpc = new G4NeutronHPCapture();
          G4NeutronHPFission* hpf = new G4NeutronHPFission();
	  Register(particle,hp,hpi,"HP");
	  Register(particle,theNeutronCapture,hpc,"HP");
	  Register(particle,theNeutronFission,hpf,"HP");
	}

        G4HadronicInteraction* theB;
        G4String s;
	if(bertFlag) {
	  theB = new G4CascadeInterface();
          s = "Bertini";
        } else if(chipsFlag) {
	  theB = new G4StringChipsInterface();
          s = "CHIPS";
	} else {
	  theB = new G4BinaryCascade();
          s = "Binary";
	}
	theB->SetMinEnergy(emin);
	theB->SetMaxEnergy(maxEcascade);
	Register(particle,hp,theB,s);
	if(glFlag) 
	  hp->AddDataSet(new G4BGGNucleonInelasticXS(particle));
	
	G4HadronicInteraction* theC = new G4LCapture();
	theC->SetMinEnergy(emin);
	Register(particle,theNeutronCapture,theC,"LCapture");

	G4HadronicInteraction* theF = new G4LFission();
	theF->SetMinEnergy(emin);
	Register(particle,theNeutronFission,theF,"LFission");

      } else if(pname == "pi-" || pname == "pi+") {
	hp->AddDataSet(&thePiCross);
        if(ftfFlag) Register(particle,hp,theFTFModel,"FTF");
	else {
	  Register(particle,hp,theQGSModel,"QGS");
	  Register(particle,hp,theFTFModel,"FTF");
	}

        if(chipsFlag) Register(particle,hp,theCHIPS,"CHIPS");
	else Register(particle,hp,theBERT,"Bertini");

	if(glFlag) 
	  hp->AddDataSet(new G4BGGPionInelasticXS(particle));

      } else if(pname == "kaon-"     || 
		pname == "kaon+"     || 
		pname == "kaon0S"    || 
		pname == "kaon0L") {
        if(ftfFlag) Register(particle,hp,theFTFModel,"FTF");
	else {
	  Register(particle,hp,theQGSModel,"QGS");
	  Register(particle,hp,theFTFModel,"FTF");
	}

        if(chipsFlag) Register(particle,hp,theCHIPS,"CHIPS");
	else Register(particle,hp,theBERT,"Bertini");

	if(glFlag) 
	  hp->AddDataSet(new G4UInelasticCrossSection(particle));

      } else if(pname == "lambda"    || 
		pname == "sigma-"    || 
		pname == "sigma+"    || 
		pname == "xi-"       || 
		pname == "xi0") {
	Register(particle,hp,theFTFModel1,"FTF");

        if(chipsFlag) Register(particle,hp,theCHIPS,"CHIPS");
	else          Register(particle,hp,theBERT,"Bertini");
	if(glFlag) 
	  hp->AddDataSet(new G4UInelasticCrossSection(particle));

      } else if(pname == "anti_proton" || pname == "anti_neutron") {
	Register(particle,hp,theFTFModel1,"FTF");
        Register(particle,hp,theCHIPS,"CHIPS");
	if(glFlag) 
	  hp->AddDataSet(new G4UInelasticCrossSection(particle));

      } else {
	Register(particle,hp,theFTFModel1,"FTF");
        Register(particle,hp,theCHIPS,"CHIPS");
	if(glFlag) 
	  hp->AddDataSet(new G4UInelasticCrossSection(particle));
      }

      if(verbose > 1)
	G4cout << "### HadronInelasticQBBC: " << hp->GetProcessName()
	       << " added for " << pname << G4endl;
    }
  }
  store->Dump(verbose);
}

void G4HadronInelasticQBBC::Register(G4ParticleDefinition* p, 
				     G4HadronicProcess* hp, 
				     G4HadronicInteraction* hi, 
				     const G4String& m)
{
  hp->RegisterMe(hi);
  store->Register(hp,p,hi,m);
  if(verbose > 1)
    G4cout << "### QBBC: Register new model " << m 
	   << " for " << p->GetParticleName() << " and " << hp->GetProcessName()
	   << " E(GeV) " << hi->GetMinEnergy()/GeV 
	   << " - " << hi->GetMaxEnergy()/GeV << G4endl;
}
