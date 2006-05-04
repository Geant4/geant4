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
// $Id: G4HadronInelasticQLHEP.cc,v 1.1 2006-05-04 16:48:39 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronInelasticQLHEP
//
// Author: 11 April 2006 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4HadronInelasticQLHEP.hh"

#include "G4HadronInelasticProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"

#include "G4PiNuclearCrossSection.hh"

#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

#include "G4CascadeInterface.hh"
#include "G4BinaryCascade.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"

#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPFission.hh"
#include "G4NeutronHPCapture.hh"

#include "G4LEAntiLambdaInelastic.hh"
#include "G4LEAntiNeutronInelastic.hh"
#include "G4LEAntiOmegaMinusInelastic.hh"
#include "G4LEAntiProtonInelastic.hh"
#include "G4LEAntiSigmaMinusInelastic.hh"
#include "G4LEAntiSigmaPlusInelastic.hh"
#include "G4LEAntiXiMinusInelastic.hh"
#include "G4LEAntiXiZeroInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LELambdaInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEOmegaMinusInelastic.hh"
#include "G4LESigmaMinusInelastic.hh"
#include "G4LESigmaPlusInelastic.hh"
#include "G4LEXiMinusInelastic.hh"
#include "G4LEXiZeroInelastic.hh"

#include "G4HadronProcessStore.hh"

G4HadronInelasticQLHEP::G4HadronInelasticQLHEP(const G4String& name, 
    G4int ver, G4bool qgs, G4bool bert, G4bool bic, G4bool hp)
  : G4VPhysicsConstructor(name), verbose(ver), qgsFlag(qgs), 
    bertFlag(bert), bicFlag(bic), hpFlag(hp), wasActivated(false)
{
  if(verbose > 1) G4cout << "### HadronInelasticQLHEP" << G4endl;
  store = G4HadronProcessStore::Instance();
  theCascade = 0;
  theQGStringDecay = 0;
  theQGStringModel = 0;
  thePreEquilib = 0;
}

G4HadronInelasticQLHEP::~G4HadronInelasticQLHEP()
{
  if(wasActivated) {
    delete theCascade;
    delete theQGStringDecay;
    delete theQGStringModel;
    delete thePreEquilib;
  }
}

void G4HadronInelasticQLHEP::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();
}

void G4HadronInelasticQLHEP::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  if(verbose > 1) G4cout << "### HadronInelasticQLHEP Construct Process" << G4endl;

  G4double minE          = 0.0*GeV;
  G4double minEcascade   = 0.0*GeV;
  if(bicFlag || bertFlag) minE = 9.5*GeV;
  if(hpFlag) minEcascade = 19.5*MeV;
  G4double minEstring    = 12.*GeV;
  G4double maxEcascade   = 25.*GeV;

  //Bertini
  G4HadronicInteraction* theBERT = 0;
  if(bertFlag) {
    theBERT = new G4CascadeInterface();
    theBERT->SetMinEnergy(minEcascade);
    theBERT->SetMaxEnergy(9.9*GeV);
  }

  //Binari
  G4HadronicInteraction* theBIC = 0;
  if(bicFlag) {
    theBIC = new G4BinaryCascade();
    theBIC->SetMinEnergy(minEcascade);
    theBIC->SetMaxEnergy(9.9*GeV);
  }

  //QGSP
  G4TheoFSGenerator* theQGSModel = 0;
  if(qgsFlag) {
    theQGSModel = new G4TheoFSGenerator();
    theQGStringModel  = new G4QGSModel< G4QGSParticipants >;
    theQGStringDecay  = new G4ExcitedStringDecay(new G4QGSMFragmentation());
    theQGStringModel->SetFragmentationModel(theQGStringDecay);
    theCascade = new G4GeneratorPrecompoundInterface;
    thePreEquilib = new G4PreCompoundModel(new G4ExcitationHandler);
    theCascade->SetDeExcitation(thePreEquilib);
    theQGSModel->SetTransport(theCascade);
    theQGSModel->SetHighEnergyGenerator(theQGStringModel);
    theQGSModel->SetMinEnergy(minEstring);
    theQGSModel->SetMaxEnergy(100*TeV);
  }

  theParticleIterator->reset();
  while( (*theParticleIterator)() ) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String pname = particle->GetParticleName();
    if(verbose > 1) G4cout << "### HadronInelasticQLHEP:  " << pname << G4endl;
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
        if(qgsFlag) {
	  Register(particle,hp,theQGSModel,"QGS");
	  AddLHEP(particle, hp, minE, maxEcascade);
	} else {
	  AddLHEP(particle, hp, minE, 100.*TeV);
	}
        if(bicFlag)       Register(particle,hp,theBIC,"Binary");
	else if(bertFlag) Register(particle,hp,theBERT,"Bertini");

      } else if(pname == "neutron") {
	hp->AddDataSet(&theXSecN);
        if(qgsFlag) {
	  Register(particle,hp,theQGSModel,"QGS");
	  AddLHEP(particle, hp, minE, maxEcascade);
	} else {
	  AddLHEP(particle, hp, minE, 100.*TeV);
	}

	G4HadronCaptureProcess* theNeutronCapture = 
	  new G4HadronCaptureProcess("nCapture");
	G4HadronFissionProcess* theNeutronFission = 
	  new G4HadronFissionProcess("nFission");
	pmanager->AddDiscreteProcess(theNeutronCapture);
	pmanager->AddDiscreteProcess(theNeutronFission);

	if(hpFlag) {
	  hp->AddDataSet(&theHPXSecI);
	  theNeutronCapture->AddDataSet(&theHPXSecC);
	  theNeutronFission->AddDataSet(&theHPXSecF);
          G4NeutronHPInelastic* hpi = new G4NeutronHPInelastic();
          G4NeutronHPCapture* hpc = new G4NeutronHPCapture();
          G4NeutronHPFission* hpf = new G4NeutronHPFission();
	  Register(particle,hp,hpi,"HP");
	  Register(particle,theNeutronCapture,hpc,"HP");
	  Register(particle,theNeutronFission,hpf,"HP");
	}

	G4HadronicInteraction* theC = new G4LCapture();
	theC->SetMinEnergy(minEcascade);
	theC->SetMaxEnergy(maxEcascade);
	Register(particle,theNeutronCapture,theC,"LCapture");

	G4HadronicInteraction* theF = new G4LFission();
	theF->SetMinEnergy(minEcascade);
	theF->SetMaxEnergy(maxEcascade);
	Register(particle,theNeutronFission,theF,"LFission");

        if(bicFlag)       Register(particle,hp,theBIC,"Binary");
	else if(bertFlag) Register(particle,hp,theBERT,"Bertini");

      } else if(pname == "pi-" || pname == "pi+") {
	hp->AddDataSet(&thePiCross);
        if(qgsFlag) {
	  Register(particle,hp,theQGSModel,"QGS");
	  AddLHEP(particle, hp, minE, maxEcascade);
	} else {
	  AddLHEP(particle, hp, minE, 100.*TeV);
	}
	if(bertFlag) Register(particle,hp,theBERT,"Bertini");

      } else if(pname == "kaon-"     || 
		pname == "kaon+"     || 
		pname == "kaon0S"    || 
		pname == "kaon0L") {

        if(qgsFlag) {
	  Register(particle,hp,theQGSModel,"QGS");
	  AddLHEP(particle, hp, minE, maxEcascade);
	} else {
	  AddLHEP(particle, hp, minE, 100.*TeV);
	}
	if(bertFlag) Register(particle,hp,theBERT,"Bertini");

      } else {

	AddLHEP(particle, hp, 0.0, 100.*TeV);
      }

      if(verbose > 1)
	G4cout << "### HadronInelasticQLHEP: " << hp->GetProcessName()
	       << " added for " << pname << G4endl;
    }
  }
}

void G4HadronInelasticQLHEP::AddLHEP(G4ParticleDefinition* particle,
				     G4HadronicProcess* hp,
				     G4double emin,
				     G4double emax)
{
  G4HadronicInteraction* hi = 0;
  G4String pname = particle->GetParticleName();

  if(pname == "anti_lambda"      ) hi = new G4LEAntiLambdaInelastic();
  else if(pname == "anti_neutron") hi = new G4LEAntiNeutronInelastic();
  else if(pname == "anti_omega-" ) hi = new G4LEAntiOmegaMinusInelastic();
  else if(pname == "anti_proton" ) hi = new G4LEAntiProtonInelastic();
  else if(pname == "anti_sigma-" ) hi = new G4LEAntiSigmaMinusInelastic();
  else if(pname == "anti_sigma+" ) hi = new G4LEAntiSigmaPlusInelastic();
  else if(pname == "anti_xi-"    ) hi = new G4LEAntiXiMinusInelastic();
  else if(pname == "anti_xi0"    ) hi = new G4LEAntiXiZeroInelastic();
  else if(pname == "kaon-"       ) hi = new G4LEKaonMinusInelastic();
  else if(pname == "kaon+"       ) hi = new G4LEKaonPlusInelastic();
  else if(pname == "kaon0S"      ) hi = new G4LEKaonZeroSInelastic();
  else if(pname == "kaon0L"      ) hi = new G4LEKaonZeroLInelastic();
  else if(pname == "neutron"     ) hi = new G4LENeutronInelastic();
  else if(pname == "lambda"      ) hi = new G4LELambdaInelastic();
  else if(pname == "omega-"      ) hi = new G4LEOmegaMinusInelastic();
  else if(pname == "proton"      ) hi = new G4LEProtonInelastic();
  else if(pname == "pi+"         ) hi = new G4LEPionPlusInelastic();
  else if(pname == "pi-"         ) hi = new G4LEPionMinusInelastic();
  else if(pname == "sigma-"      ) hi = new G4LESigmaMinusInelastic();
  else if(pname == "sigma+"      ) hi = new G4LESigmaPlusInelastic();
  else if(pname == "xi-"         ) hi = new G4LEXiMinusInelastic();
  else if(pname == "xi0"         ) hi = new G4LEXiZeroInelastic();

  if(hi) {
    hp->RegisterMe(hi);
    hi->SetMinEnergy(emin);
    hi->SetMaxEnergy(emax);
    Register(particle,hp,hi,"LHEP");
  } else {
    G4cout << "### G4HadronInelasticTHEO: ERROR - no LHEP model for "
           << pname << G4endl;
  }
}

void G4HadronInelasticQLHEP::Register(G4ParticleDefinition* p, G4HadronicProcess* hp, 
				      G4HadronicInteraction* hi, const G4String& m)
{
  hp->RegisterMe(hi);
  store->Register(hp,p,hi,m);
  if(verbose > 1)
    G4cout << "### QLHEP: Register new model " << m 
	   << " for " << p->GetParticleName() << " and " << hp->GetProcessName()
	   << " E(GeV) " << hi->GetMinEnergy()/GeV 
	   << " - " << hi->GetMaxEnergy()/GeV << G4endl;
}
