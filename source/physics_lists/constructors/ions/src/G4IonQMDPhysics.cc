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
// $Id: G4IonQMDPhysics.cc 80671 2014-05-06 13:59:16Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4IonQMDPhysics
//   Created from G4IonBinaryCascadePhysics
//
// Author:      G.Folger
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4IonQMDPhysics.hh"

#include "G4SystemOfUnits.hh"

#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"
#include "G4IonConstructor.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4QMDReaction.hh"

#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4FTFBuilder.hh"
#include "G4HadronicInteraction.hh"
#include "G4BuilderType.hh"

#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionInelastic.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

// Nuclei
#include "G4IonConstructor.hh"
#include "G4BuilderType.hh"
#include "G4HadronicInteractionRegistry.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4IonQMDPhysics);

G4ThreadLocal std::vector<G4HadronInelasticProcess*>* G4IonQMDPhysics::p_list = 0;
G4ThreadLocal std::vector<G4HadronicInteraction*>* G4IonQMDPhysics::model_list = 0;

G4ThreadLocal G4VCrossSectionDataSet* G4IonQMDPhysics::theNuclNuclData = 0; 
G4ThreadLocal G4VComponentCrossSection* G4IonQMDPhysics::theGGNuclNuclXS = 0;

G4ThreadLocal G4BinaryLightIonReaction* G4IonQMDPhysics::theIonBC = 0;
G4ThreadLocal G4HadronicInteraction* G4IonQMDPhysics::theFTFP = 0;
G4ThreadLocal G4FTFBuilder* G4IonQMDPhysics::theBuilder = 0;
G4ThreadLocal  G4QMDReaction* G4IonQMDPhysics::theQMD = 0;
G4ThreadLocal G4bool G4IonQMDPhysics::wasActivated = false;

G4IonQMDPhysics::G4IonQMDPhysics(G4int ver)
  :  G4VPhysicsConstructor("IonQMD"), verbose(ver)
{
  eminBIC  = 0.*MeV;
  eminQMD  = 100.*MeV;
  emaxQMD  = 10.*GeV;
  emaxFTFP = 1.*TeV;
  overlap  = 10*MeV;
  SetPhysicsType(bIons);
  if(verbose > 1) G4cout << "### G4IonQMDPhysics" << G4endl;
}

G4IonQMDPhysics::G4IonQMDPhysics(const G4String& name, 
						     G4int ver)
  :  G4VPhysicsConstructor(name), verbose(ver)
{
  eminBIC  = 0.*MeV;
  eminQMD  = 100.*MeV;
  emaxQMD  = 10.*GeV;
  emaxFTFP = 1.*TeV;
  overlap  = 10*MeV;
  SetPhysicsType(bIons);
  if(verbose > 1) G4cout << "### G4IonQMDPhysics" << G4endl;
}

G4IonQMDPhysics::~G4IonQMDPhysics()
{
  if(wasActivated) {
    delete theBuilder; theBuilder = 0;
    delete theGGNuclNuclXS; theGGNuclNuclXS = 0; 
    delete theNuclNuclData; theNuclNuclData = 0;
    G4int i;
    if ( p_list ) {
      G4int n = p_list->size();
      for(i=0; i<n; i++) {delete (*p_list)[i];}
      delete p_list; p_list = 0;
    }
    if ( model_list ) {
      G4int n = model_list->size();
      for(i=0; i<n; i++) {delete (*model_list)[i];}
      delete model_list; model_list = 0;
    }
  }
}

void G4IonQMDPhysics::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
  G4PreCompoundModel* thePreCompound = static_cast<G4PreCompoundModel*>(p); 
  if(!thePreCompound) { thePreCompound = new G4PreCompoundModel; }

  theIonBC = new G4BinaryLightIonReaction(thePreCompound);
  if ( model_list == 0 ) model_list = new std::vector<G4HadronicInteraction*>;
  model_list->push_back(theIonBC);

  theBuilder = new G4FTFBuilder("FTFP",thePreCompound);
  theFTFP = theBuilder->GetModel();
  model_list->push_back(theFTFP);

  theQMD= new G4QMDReaction();
  model_list->push_back(theQMD);

  theNuclNuclData = new G4CrossSectionInelastic( theGGNuclNuclXS = new G4ComponentGGNuclNuclXsc() );

  AddProcess("dInelastic", G4Deuteron::Deuteron(), theIonBC, theQMD, theFTFP);
  AddProcess("tInelastic", G4Triton::Triton(), theIonBC, theQMD, theFTFP);
  AddProcess("He3Inelastic", G4He3::He3(), theIonBC, theQMD, theFTFP);
  AddProcess("alphaInelastic", G4Alpha::Alpha(), theIonBC, theQMD, theFTFP);
  AddProcess("ionInelastic", G4GenericIon::GenericIon(), theIonBC, theQMD, theFTFP);

}

void G4IonQMDPhysics::AddProcess(const G4String& name,
			         G4ParticleDefinition* p, 
				 G4BinaryLightIonReaction* BIC,
				 G4QMDReaction* QMD,
				 G4HadronicInteraction* FTFP)
{
  G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess(name, p);
  if ( p_list == 0 ) p_list = new  std::vector<G4HadronInelasticProcess*>;
  p_list->push_back(hadi);
  G4ProcessManager* pManager = p->GetProcessManager();
  pManager->AddDiscreteProcess(hadi);

  hadi->AddDataSet(theNuclNuclData);

  BIC->SetMinEnergy(eminBIC);
  BIC->SetMaxEnergy(emaxQMD);  //reset when QMD is present 
  hadi->RegisterMe(BIC);

  if(QMD) {
    QMD->SetMinEnergy(eminQMD);
    BIC->SetMaxEnergy(eminQMD+overlap);
    QMD->SetMaxEnergy(emaxQMD);
    hadi->RegisterMe(QMD);
  }  

  if(FTFP) {
    FTFP->SetMinEnergy(emaxQMD - overlap);
    FTFP->SetMaxEnergy(emaxFTFP);
    hadi->RegisterMe(FTFP);
  }  

  if(verbose > 1) {
    G4cout << "Register " << hadi->GetProcessName()
	   << " for " << p->GetParticleName() << G4endl
	   << "       Binary Cascade for E(MeV)= " << eminBIC << " - " 
	   << (QMD==0 ? emaxQMD : (eminQMD-overlap)) ;
    if(QMD) {
      G4cout  << G4endl <<"       QMD for E(MeV)= " << eminQMD << " - " << emaxQMD;
    }
    if(FTFP) {
      G4cout << G4endl<< "       FTFP for E(MeV)= " << emaxQMD-overlap << " - " << emaxFTFP;
    }
    G4cout << G4endl;
  }
}

void G4IonQMDPhysics::ConstructParticle()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}
