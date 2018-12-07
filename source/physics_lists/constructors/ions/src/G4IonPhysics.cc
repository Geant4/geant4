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
// Class:    G4IonPhysics
//
// Author:      A.Ivanchenko 02.03.2011
//
// Modified: 
// 16.10.2012 A.Ribon: renamed G4IonFTFPBinaryCascadePhysics as G4IonPhysics     
//
//---------------------------------------------------------------------------

#include "G4IonPhysics.hh"
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

#include "G4HadronicParameters.hh"

using namespace std;

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4IonPhysics);

G4ThreadLocal G4FTFBuilder* G4IonPhysics::theBuilder=nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4IonPhysics::G4IonPhysics(G4int ver)
  : G4IonPhysics("ionInelasticFTFP_BIC")
{
  verbose = ver;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4IonPhysics::G4IonPhysics(const G4String& nname)
  : G4VPhysicsConstructor(nname),verbose(1)
{
  SetPhysicsType(bIons);
  if(verbose > 1) { G4cout << "### IonPhysics: " << nname << G4endl; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4IonPhysics::~G4IonPhysics()
{
  // Explictly setting pointers to zero is actually needed.
  // These are static variables, in case we restart threads 
  // we need to re-create objects
  delete theBuilder; 
  theBuilder = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4IonPhysics::ConstructParticle()
{
  //  Construct ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4IonPhysics::ConstructProcess()
{
  G4double emax = G4HadronicParameters::Instance()->GetMaxEnergy();

  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
  G4PreCompoundModel* thePreCompound = static_cast<G4PreCompoundModel*>(p); 
  if(!thePreCompound) { thePreCompound = new G4PreCompoundModel; }

  // Binary Cascade
  G4HadronicInteraction* theIonBC = 
    new G4BinaryLightIonReaction(thePreCompound);
  theIonBC->SetMinEnergy(0.0);
  theIonBC->SetMaxEnergy(4*GeV);

  // FTFP
  G4HadronicInteraction* theFTFP = nullptr;
  if(emax > theIonBC->GetMaxEnergy()) {
    theBuilder = new G4FTFBuilder("FTFP",thePreCompound);
    theFTFP = theBuilder->GetModel();
    theFTFP->SetMinEnergy(2*GeV);
    theFTFP->SetMaxEnergy(emax);
  }

  G4CrossSectionInelastic* theNuclNuclData = 
    new G4CrossSectionInelastic(new G4ComponentGGNuclNuclXsc());

  AddProcess("dInelastic", G4Deuteron::Deuteron(),theIonBC,theFTFP,
	     theNuclNuclData);
  AddProcess("tInelastic",G4Triton::Triton(),theIonBC,theFTFP,
	     theNuclNuclData);
  AddProcess("He3Inelastic",G4He3::He3(),theIonBC,theFTFP,
	     theNuclNuclData);
  AddProcess("alphaInelastic", G4Alpha::Alpha(),theIonBC,theFTFP,
	     theNuclNuclData);
  AddProcess("ionInelastic",G4GenericIon::GenericIon(),theIonBC,theFTFP,
	     theNuclNuclData);

  if(verbose > 1) {
    G4cout << "G4IonPhysics::ConstructProcess done! " 
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4IonPhysics::AddProcess(const G4String& name, 
			      G4ParticleDefinition* part,
			      G4HadronicInteraction* theIonBC,
			      G4HadronicInteraction* theFTFP,
			      G4VCrossSectionDataSet* xs)
{
  G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess(name, part);
  G4ProcessManager* pManager = part->GetProcessManager();
  pManager->AddDiscreteProcess(hadi);
    
  hadi->AddDataSet(xs);
    
  hadi->RegisterMe(theIonBC);
  if(theFTFP) { hadi->RegisterMe(theFTFP); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
