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
// $Id: G4HadronInelasticQLHEP.cc,v 1.4 2010-06-08 08:58:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronInelasticQLHEP
//
// Author: 11 April 2006 V. Ivanchenko
//
// Modified:
// 05.07.2006 V.Ivanchenko fix problem of initialisation of HP
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

#include "G4HEAntiLambdaInelastic.hh"
#include "G4HEAntiNeutronInelastic.hh"
#include "G4HEAntiOmegaMinusInelastic.hh"
#include "G4HEAntiProtonInelastic.hh"
#include "G4HEAntiSigmaMinusInelastic.hh"
#include "G4HEAntiSigmaPlusInelastic.hh"
#include "G4HEAntiXiMinusInelastic.hh"
#include "G4HEAntiXiZeroInelastic.hh"
#include "G4HEKaonMinusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HENeutronInelastic.hh"
#include "G4HELambdaInelastic.hh"
#include "G4HEProtonInelastic.hh"
#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4HEOmegaMinusInelastic.hh"
#include "G4HESigmaMinusInelastic.hh"
#include "G4HESigmaPlusInelastic.hh"
#include "G4HEXiMinusInelastic.hh"
#include "G4HEXiZeroInelastic.hh"

G4HadronInelasticQLHEP::G4HadronInelasticQLHEP(G4int ver)
  : G4VPhysicsConstructor("hInelasticQLHEP"), verbose(ver), qgsFlag(false), 
    bertFlag(false), bicFlag(false), hpFlag(false), wasActivated(false)
{
  if(verbose > 1) G4cout << "### HadronInelasticQLHEP" << G4endl;
  theCascade = 0;
  theQGStringDecay = 0;
  theQGStringModel = 0;
  thePreEquilib = 0;
  theHPXSecI = 0;
  theHPXSecC = 0;
  theHPXSecF = 0;
}

G4HadronInelasticQLHEP::G4HadronInelasticQLHEP(const G4String&, 
    G4int ver, G4bool qgs, G4bool bert, G4bool bic, G4bool hp)
  : G4VPhysicsConstructor("hInelasticQLHEP"), verbose(ver), qgsFlag(qgs), 
    bertFlag(bert), bicFlag(bic), hpFlag(hp), wasActivated(false)
{
  if(verbose > 1) G4cout << "### HadronInelasticQLHEP" << G4endl;
  theCascade = 0;
  theQGStringDecay = 0;
  theQGStringModel = 0;
  thePreEquilib = 0;
  theHPXSecI = 0;
  theHPXSecC = 0;
  theHPXSecF = 0;
}

G4HadronInelasticQLHEP::~G4HadronInelasticQLHEP()
{
  if(wasActivated) {
    delete theCascade;
    delete theQGStringDecay;
    delete theQGStringModel;
    delete thePreEquilib;
    delete theHPXSecI;
    delete theHPXSecC;
    delete theHPXSecF;
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

  G4double minEneutron   = 0.0*GeV;
  G4double minELEP       = 0.0*GeV;
  G4double maxELEP       = 55.*GeV;
  G4double minEstring    = 12.*GeV;
  if(hpFlag) minEneutron = 19.5*MeV;

  //Bertini
  G4HadronicInteraction* theBERT = 0;
  if(bertFlag) {
    minELEP = 9.5*GeV;
    theBERT = new G4CascadeInterface();
    theBERT->SetMinEnergy(0.0);
    theBERT->SetMaxEnergy(9.9*GeV);
  }

  //Binari
  G4HadronicInteraction* theBIC = 0;
  if(bicFlag) {
    minELEP = 9.5*GeV;
    theBIC = new G4BinaryCascade();
    theBIC->SetMinEnergy(0.0);
    theBIC->SetMaxEnergy(9.9*GeV);
  }

  //QGSP
  G4TheoFSGenerator* theQGSModel = 0;
  if(qgsFlag) {
    maxELEP     = 25.*GeV;
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
        if(qgsFlag) {
	  Register(particle,hp,theQGSModel,"QGS");
	  hp->AddDataSet(&theXSecP);
	} else {
	  AddHEP(particle, hp, 25.*GeV, 100.*TeV);
	}
	AddLEP(particle, hp, minELEP, maxELEP);
       
        if(bicFlag)       Register(particle,hp,theBIC,"Binary");
	else if(bertFlag) Register(particle,hp,theBERT,"Bertini");

      } else if(pname == "neutron") {
        if(qgsFlag) {
	  Register(particle,hp,theQGSModel,"QGS");
	  hp->AddDataSet(&theXSecN);
	} else {
	  AddHEP(particle, hp, 25.*GeV, 100.*TeV);
	}
	AddLEP(particle, hp, minELEP, maxELEP);

	G4HadronCaptureProcess* theNeutronCapture = 
	  new G4HadronCaptureProcess("nCapture");
	G4HadronFissionProcess* theNeutronFission = 
	  new G4HadronFissionProcess("nFission");
	pmanager->AddDiscreteProcess(theNeutronCapture);
	pmanager->AddDiscreteProcess(theNeutronFission);

	if(hpFlag) {
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

	G4HadronicInteraction* theC = new G4LCapture();
	theC->SetMinEnergy(minEneutron);
	Register(particle,theNeutronCapture,theC,"LCapture");

	G4HadronicInteraction* theF = new G4LFission();
	theF->SetMinEnergy(minEneutron);
	Register(particle,theNeutronFission,theF,"LFission");

        if(bicFlag) {
	  G4BinaryCascade* theB = new G4BinaryCascade();
	  theB->SetMinEnergy(minEneutron);
	  theB->SetMaxEnergy(9.9*GeV);
	  Register(particle,hp,theB,"Binary");
	} else if(bertFlag) {
	  G4CascadeInterface* theB = new G4CascadeInterface();
	  theB->SetMinEnergy(minEneutron);
	  theB->SetMaxEnergy(9.9*GeV);
	  Register(particle,hp,theB,"Bertini");
	}

      } else if(pname == "pi-" || pname == "pi+") {
        if(qgsFlag) {
	  Register(particle,hp,theQGSModel,"QGS");
	  hp->AddDataSet(&thePiCross);
	} else {
	  AddHEP(particle, hp, 25.*GeV, 100.*TeV);
	}
	AddLEP(particle, hp, minELEP, maxELEP);
	if(bertFlag) Register(particle,hp,theBERT,"Bertini");

      } else if(pname == "kaon-"     || 
		pname == "kaon+"     || 
		pname == "kaon0S"    || 
		pname == "kaon0L") {

        if(qgsFlag) Register(particle,hp,theQGSModel,"QGS");
	else        AddHEP(particle, hp, 25.*GeV, 100.*TeV);
       
	AddLEP(particle, hp, minELEP, maxELEP);
	if(bertFlag) Register(particle,hp,theBERT,"Bertini");

      } else {

	AddHEP(particle, hp, 25.*GeV, 100.*TeV);
	AddLEP(particle, hp, 0.0, 55.*GeV);
      }

      if(verbose > 1)
	G4cout << "### HadronInelasticQLHEP: " << hp->GetProcessName()
	       << " added for " << pname << G4endl;
    }
  }
}

void G4HadronInelasticQLHEP::AddLEP(G4ParticleDefinition* particle,
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
    hi->SetMinEnergy(emin);
    hi->SetMaxEnergy(emax);
    Register(particle,hp,hi,"LHEP");
  } else {
    G4cout << "### G4HadronInelasticTHEO: ERROR - no LHEP model for "
           << pname << G4endl;
  }
}


void G4HadronInelasticQLHEP::AddHEP(G4ParticleDefinition* particle,
				    G4HadronicProcess* hp,
				    G4double emin,
				    G4double emax)
{
  G4HadronicInteraction* hi = 0;
  G4String pname = particle->GetParticleName();

  if(pname == "anti_lambda"      ) hi = new G4HEAntiLambdaInelastic();
  else if(pname == "anti_neutron") hi = new G4HEAntiNeutronInelastic();
  else if(pname == "anti_omega-" ) hi = new G4HEAntiOmegaMinusInelastic();
  else if(pname == "anti_proton" ) hi = new G4HEAntiProtonInelastic();
  else if(pname == "anti_sigma-" ) hi = new G4HEAntiSigmaMinusInelastic();
  else if(pname == "anti_sigma+" ) hi = new G4HEAntiSigmaPlusInelastic();
  else if(pname == "anti_xi-"    ) hi = new G4HEAntiXiMinusInelastic();
  else if(pname == "anti_xi0"    ) hi = new G4HEAntiXiZeroInelastic();
  else if(pname == "kaon-"       ) hi = new G4HEKaonMinusInelastic();
  else if(pname == "kaon+"       ) hi = new G4HEKaonPlusInelastic();
  else if(pname == "kaon0S"      ) hi = new G4HEKaonZeroInelastic();
  else if(pname == "kaon0L"      ) hi = new G4HEKaonZeroInelastic();
  else if(pname == "neutron"     ) hi = new G4HENeutronInelastic();
  else if(pname == "lambda"      ) hi = new G4HELambdaInelastic();
  else if(pname == "omega-"      ) hi = new G4HEOmegaMinusInelastic();
  else if(pname == "proton"      ) hi = new G4HEProtonInelastic();
  else if(pname == "pi+"         ) hi = new G4HEPionPlusInelastic();
  else if(pname == "pi-"         ) hi = new G4HEPionMinusInelastic();
  else if(pname == "sigma-"      ) hi = new G4HESigmaMinusInelastic();
  else if(pname == "sigma+"      ) hi = new G4HESigmaPlusInelastic();
  else if(pname == "xi-"         ) hi = new G4HEXiMinusInelastic();
  else if(pname == "xi0"         ) hi = new G4HEXiZeroInelastic();

  if(hi) {
    hi->SetMinEnergy(emin);
    hi->SetMaxEnergy(emax);
    Register(particle,hp,hi,"HEP");
  } else {
    G4cout << "### G4HadronInelasticQLHEP: ERROR - no HEP model for "
           << pname << G4endl;
  }
}

void G4HadronInelasticQLHEP::Register(G4ParticleDefinition* p, 
				      G4HadronicProcess* hp, 
				      G4HadronicInteraction* hi, 
				      const G4String& m)
{
  hp->RegisterMe(hi);
  if(verbose > 1)
    G4cout << "### QLHEP: Register new model " << m 
	   << " for " << p->GetParticleName() << " and " 
	   << hp->GetProcessName()
	   << " E(GeV) " << hi->GetMinEnergy()/GeV 
	   << " - " << hi->GetMaxEnergy()/GeV << G4endl;
}
