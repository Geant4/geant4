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
// $Id: G4IonBinaryCascadePhysics.cc 80671 2014-05-06 13:59:16Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4IonBinaryCascadePhysics
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
// 21.03.13 A.Ribon : replace LHEP with FTFP.
// 23.06.06 V.Ivanchenko set emaxLHEP=1 TeV
// 24.06.06 V.Ivanchenko fix typo
//
//----------------------------------------------------------------------------
//

#include "G4IonBinaryCascadePhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"
#include "G4IonConstructor.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionInelastic.hh"

#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4FTFBuilder.hh"
#include "G4HadronicInteraction.hh"
#include "G4BuilderType.hh"
#include "G4HadronicInteractionRegistry.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4IonBinaryCascadePhysics);

G4ThreadLocal G4bool G4IonBinaryCascadePhysics::wasActivated = false;
G4ThreadLocal std::vector<G4HadronInelasticProcess*>* G4IonBinaryCascadePhysics::G4MT_p_list = 0;
G4ThreadLocal std::vector<G4HadronicInteraction*>* G4IonBinaryCascadePhysics::G4MT_model_list = 0;

G4ThreadLocal G4VCrossSectionDataSet* G4IonBinaryCascadePhysics::theNuclNuclData = 0; 
G4ThreadLocal G4VComponentCrossSection* G4IonBinaryCascadePhysics::theGGNuclNuclXS = 0;
G4ThreadLocal G4BinaryLightIonReaction* G4IonBinaryCascadePhysics::theIonBC = 0;
G4ThreadLocal G4HadronicInteraction* G4IonBinaryCascadePhysics::theFTFP = 0;
G4ThreadLocal G4FTFBuilder* G4IonBinaryCascadePhysics::theBuilder = 0;

G4IonBinaryCascadePhysics::G4IonBinaryCascadePhysics(G4int ver)
  :  G4VPhysicsConstructor("IonBinaryCascade"), verbose(ver)
{
  theNuclNuclData = 0; 
  theGGNuclNuclXS = 0;
  theIonBC = 0;
  theFTFP = 0;
  theBuilder = 0;
  SetPhysicsType(bIons);
  if(verbose > 1) G4cout << "### G4IonBinaryCascadePhysics" << G4endl;
}


G4IonBinaryCascadePhysics::G4IonBinaryCascadePhysics(const G4String& name, 
						     G4int ver)
  :  G4VPhysicsConstructor(name), verbose(ver)
{
  theNuclNuclData = 0; 
  theGGNuclNuclXS = 0;
  theIonBC = 0;
  theFTFP = 0;
  theBuilder = 0;
  SetPhysicsType(bIons);
  if(verbose > 1) G4cout << "### G4IonBinaryCascadePhysics" << G4endl;
}


G4IonBinaryCascadePhysics::~G4IonBinaryCascadePhysics()
{
  //Explicitly setting pointers TLS to 0 is 
  //needed, in case I create a new thread
  //this variable is static!
  if(wasActivated) {
    delete theBuilder; theBuilder = 0;
    delete theGGNuclNuclXS; theGGNuclNuclXS = 0;
    delete theNuclNuclData;  theNuclNuclData = 0;
    G4int i;
    if ( G4MT_p_list ) {
      G4int n = G4MT_p_list->size();
      for(i=0; i<n; i++) {delete (*G4MT_p_list)[i];}
      delete G4MT_p_list;
      G4MT_p_list = 0; 
    }
    if ( G4MT_model_list ) {
      G4int n = G4MT_model_list->size();
      for(i=0; i<n; i++) {delete (*G4MT_model_list)[i];}
      delete G4MT_model_list;
      G4MT_model_list = 0;
    }
  }
}


void G4IonBinaryCascadePhysics::ConstructProcess()
{
  if(wasActivated) { return; }
  wasActivated = true;
  
  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
  G4PreCompoundModel* thePreCompound = static_cast<G4PreCompoundModel*>(p); 
  if(!thePreCompound) { thePreCompound = new G4PreCompoundModel; }

  theIonBC = new G4BinaryLightIonReaction(thePreCompound);
  theIonBC->SetMinEnergy(0.0);
  theIonBC->SetMaxEnergy(4.0*GeV);
  if ( G4MT_model_list == 0 ) G4MT_model_list = new std::vector<G4HadronicInteraction*>;
  G4MT_model_list->push_back(theIonBC);

  theBuilder = new G4FTFBuilder("FTFP",thePreCompound);
  theFTFP = theBuilder->GetModel();
  theFTFP->SetMinEnergy(2.0*GeV);
  theFTFP->SetMaxEnergy(100.0*TeV);
  G4MT_model_list->push_back(theFTFP);

  theNuclNuclData = new G4CrossSectionInelastic( theGGNuclNuclXS = new G4ComponentGGNuclNuclXsc() );

  AddProcess("dInelastic", G4Deuteron::Deuteron());
  AddProcess("tInelastic", G4Triton::Triton());
  AddProcess("He3Inelastic", G4He3::He3());
  AddProcess("alphaInelastic", G4Alpha::Alpha());
  AddProcess("ionInelastic", G4GenericIon::GenericIon());
}


void G4IonBinaryCascadePhysics::AddProcess(const G4String& name, 
                                           G4ParticleDefinition* part)
{
  if ( G4MT_p_list == 0 ) G4MT_p_list = new std::vector<G4HadronInelasticProcess*>;
  G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess(name, part);
  G4MT_p_list->push_back(hadi);
  G4ProcessManager* pManager = part->GetProcessManager();
  pManager->AddDiscreteProcess(hadi);

  hadi->AddDataSet(theNuclNuclData);
    
  hadi->RegisterMe(theIonBC);
  hadi->RegisterMe(theFTFP);
}


void G4IonBinaryCascadePhysics::ConstructParticle()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}
