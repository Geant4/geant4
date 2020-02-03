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
//
//---------------------------------------------------------------------------
//
// ClassName:   G4IonINCLXXPhysics
//
// Author:      D. Mancusi
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4IonINCLXXPhysics.hh"

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
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionInelastic.hh"

#include "G4INCLXXInterface.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4FTFBuilder.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4HadronicParameters.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclearLevelData.hh"
#include "G4BuilderType.hh"
#include "G4HadronicParameters.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4IonINCLXXPhysics);

G4ThreadLocal G4INCLXXInterface* G4IonINCLXXPhysics::theINCLXXDeuteron = nullptr;
G4ThreadLocal G4INCLXXInterface* G4IonINCLXXPhysics::theINCLXXTriton = nullptr;
G4ThreadLocal G4INCLXXInterface* G4IonINCLXXPhysics::theINCLXXHe3 = nullptr;
G4ThreadLocal G4INCLXXInterface* G4IonINCLXXPhysics::theINCLXXAlpha = nullptr;
G4ThreadLocal G4INCLXXInterface* G4IonINCLXXPhysics::theINCLXXIons = nullptr;
G4ThreadLocal G4FTFBuilder* G4IonINCLXXPhysics::theFTFPBuilder = nullptr;

G4IonINCLXXPhysics::G4IonINCLXXPhysics(G4int ver) :
  G4IonINCLXXPhysics("IonINCLXX", ver)
{}

G4IonINCLXXPhysics::G4IonINCLXXPhysics(const G4String& nname, G4int ver)
  : G4VPhysicsConstructor(nname), verbose(ver)
{
  // INCLXX light ion maximum energy is 3.0 GeV/nucleon
  emaxINCLXX = 3.0 * GeV;
  deltaE     = 100.*MeV;
  SetPhysicsType(bIons);
  G4DeexPrecoParameters* param = G4NuclearLevelData::GetInstance()->GetParameters();
  param->SetDeexChannelsType(fCombined);
  if(verbose > 1) { G4cout << "### IonPhysics: " << nname << G4endl; }
}

G4IonINCLXXPhysics::~G4IonINCLXXPhysics()
{
  //For MT need to explicitly set back pointers to zero:
  //variables are static and if new threads are created we can have problems
  //since variable is still pointing old value
  delete theFTFPBuilder; theFTFPBuilder=nullptr;
  delete theINCLXXDeuteron; theINCLXXDeuteron = nullptr;
  delete theINCLXXTriton; theINCLXXTriton=nullptr;
  delete theINCLXXHe3; theINCLXXHe3=nullptr;
  delete theINCLXXAlpha; theINCLXXAlpha=nullptr;
  delete theINCLXXIons; theINCLXXIons=nullptr;
}

void G4IonINCLXXPhysics::ConstructProcess()
{
  theINCLXXDeuteron= new G4INCLXXInterface();
  theINCLXXTriton= new G4INCLXXInterface();
  theINCLXXHe3= new G4INCLXXInterface();
  theINCLXXAlpha= new G4INCLXXInterface();
  theINCLXXIons= new G4INCLXXInterface();

  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
  G4PreCompoundModel* thePreCompound = static_cast<G4PreCompoundModel*>(p); 
  if(!thePreCompound) { thePreCompound = new G4PreCompoundModel; }

  G4CrossSectionInelastic* theNuclNuclData = 
    new G4CrossSectionInelastic( new G4ComponentGGNuclNuclXsc() );

  G4double emax = G4HadronicParameters::Instance()->GetMaxEnergy();
  G4HadronicInteraction* theFTFP = nullptr;
  if(emax > emaxINCLXX) {
    theFTFPBuilder = new G4FTFBuilder("FTFP",thePreCompound);
    theFTFP = theFTFPBuilder->GetModel();
    theFTFP->SetMinEnergy(emaxINCLXX - deltaE);
    theFTFP->SetMaxEnergy(emax);
  }

  AddProcess("dInelastic", G4Deuteron::Deuteron(), theINCLXXDeuteron, theFTFP, theNuclNuclData);
  AddProcess("tInelastic", G4Triton::Triton(), theINCLXXTriton, theFTFP, theNuclNuclData);
  AddProcess("He3Inelastic", G4He3::He3(), theINCLXXHe3, theFTFP, theNuclNuclData);
  AddProcess("alphaInelastic", G4Alpha::Alpha(), theINCLXXAlpha, theFTFP, theNuclNuclData);
  AddProcess("ionInelastic", G4GenericIon::GenericIon(), theINCLXXIons, theFTFP, theNuclNuclData);
}

void G4IonINCLXXPhysics::AddProcess(const G4String& name,
				    G4ParticleDefinition* p, 
				    G4HadronicInteraction* model1,
				    G4HadronicInteraction* model2,
				    G4VCrossSectionDataSet* xs)
{
  G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess(name, p);
  G4ProcessManager* pManager = p->GetProcessManager();
  pManager->AddDiscreteProcess(hadi);
  hadi->AddDataSet(xs);    
  model1->SetMaxEnergy(emaxINCLXX);
  hadi->RegisterMe(model1);
  if(model2) { hadi->RegisterMe(model2); }

  if(verbose > 1) {
    G4cout << "Register " << hadi->GetProcessName()
	   << " for " << p->GetParticleName()
	   << " INCLXX/G4DeexcitationHandler for E(MeV)= 0" << " - " << emaxINCLXX;
    if(model2) {
      G4cout  << " FTFP for E(MeV)= " << emaxINCLXX - deltaE << " - " 
	      << model2->GetMaxEnergy();
    }
    G4cout << G4endl;
  }
}

void G4IonINCLXXPhysics::ConstructParticle()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}
