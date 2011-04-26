// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software.                                *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This code implementation is the intellectual property of the ESA.*
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4IonFTFPBinaryCascadePhysics.cc,v 1.1 2006/10/28 16:00:25 vnivanch Exp $
// GEANT4 tag $Name: $
//
//---------------------------------------------------------------------------
//
// Class:    G4IonFTFPBinaryCascadePhysics
//
// Author:      A.Ivanchenko 02.03.2011
//
//
// Modified:
//
// ------------------------------------------------------------
// 

#include "G4IonFTFPBinaryCascadePhysics.hh"
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
#include "G4TripathiCrossSection.hh"
#include "G4TripathiLightCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4IonProtonCrossSection.hh"

#include "G4PreCompoundModel.hh"
#include "G4FTFBuilder.hh"
#include "G4HadronicInteraction.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4IonFTFPBinaryCascadePhysics::G4IonFTFPBinaryCascadePhysics(G4int ver)
  : G4VPhysicsConstructor("ionInelasticFTFP_BIC"),verbose(ver),
    wasActivated(false)
{
  fTripathi = fTripathiLight = fShen = fIonH = 0;
  theIonBC = 0;
  theFTFP = 0;
  theBuilder = 0;
  if(verbose > 1) { G4cout << "### G4IonFTFPBinaryCascadePhysics" << G4endl; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4IonFTFPBinaryCascadePhysics::~G4IonFTFPBinaryCascadePhysics()
{
  delete theBuilder;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4IonFTFPBinaryCascadePhysics::ConstructParticle()
{
  //  Construct ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4IonFTFPBinaryCascadePhysics::ConstructProcess()
{
  if(wasActivated) { return; }
  wasActivated = true;

  G4double emax = 100.*TeV;

  G4PreCompoundModel* thePreCompound = 
    new G4PreCompoundModel(new G4ExcitationHandler());

  // Binary Cascade
  theIonBC = new G4BinaryLightIonReaction();
  theIonBC->SetPrecompound(thePreCompound);
  theIonBC->SetMinEnergy(0.0);
  theIonBC->SetMaxEnergy(6*GeV);

  // FTFP
  theBuilder = new G4FTFBuilder("FTFP",thePreCompound);
  theFTFP = theBuilder->GetModel();
  theFTFP->SetMinEnergy(5*GeV);
  theFTFP->SetMaxEnergy(emax);

  fShen = new G4IonsShenCrossSection();
  fTripathi = new G4TripathiCrossSection();
  fTripathiLight = new G4TripathiLightCrossSection();
  fIonH = new G4IonProtonCrossSection();

  AddProcess("dInelastic", G4Deuteron::Deuteron(),false);
  AddProcess("tInelastic",G4Triton::Triton(),false);
  AddProcess("He3Inelastic",G4He3::He3(),true);
  AddProcess("alphaInelastic", G4Alpha::Alpha(),true);
  AddProcess("ionInelastic",G4GenericIon::GenericIon(),true);

  if(verbose > 1) {
    G4cout << "G4IonFTFPBinaryCascadePhysics::ConstructProcess done! " 
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4IonFTFPBinaryCascadePhysics::AddProcess(const G4String& name, 
					       G4ParticleDefinition* part, 
					       G4bool isIon)
{
  G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess(name, part);
  G4ProcessManager* pManager = part->GetProcessManager();
  pManager->AddDiscreteProcess(hadi);
  hadi->AddDataSet(fShen);
  hadi->AddDataSet(fTripathi);
  hadi->AddDataSet(fTripathiLight);
  if(isIon) { hadi->AddDataSet(fIonH); }
  hadi->RegisterMe(theIonBC);
  hadi->RegisterMe(theFTFP);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
