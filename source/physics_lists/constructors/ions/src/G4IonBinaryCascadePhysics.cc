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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4IonBinaryCascadePhysics
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
// 23.06.06 V.Ivanchenko set emaxLHEP=1 TeV
// 24.06.06 V.Ivanchenko fix typo
//
//----------------------------------------------------------------------------
//

#include "G4IonBinaryCascadePhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"

#include "G4BinaryLightIonReaction.hh"
#include "G4TripathiCrossSection.hh"
#include "G4TripathiLightCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4IonProtonCrossSection.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

// Nuclei
#include "G4IonConstructor.hh"
#include "G4BuilderType.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4IonBinaryCascadePhysics);


G4IonBinaryCascadePhysics::G4IonBinaryCascadePhysics(G4int ver)
  :  G4VPhysicsConstructor("IonBinaryCascade"), verbose(ver), wasActivated(false)
{
  fLEDModel = 0;
  fLETModel = 0;
  fLEAModel = 0;
  fTripathi = 0; 
  fTripathiLight = 0;
  fShen = 0;
  fIonH = 0;
  emax     = 20.*GeV;
  emaxLHEP = 1.*TeV;
  eminBIC  = 0.*MeV;
  SetPhysicsType(bIons);
  if(verbose > 1) G4cout << "### G4IonBinaryCascadePhysics" << G4endl;
}

G4IonBinaryCascadePhysics::G4IonBinaryCascadePhysics(const G4String& name, 
						     G4int ver)
  :  G4VPhysicsConstructor(name), verbose(ver), wasActivated(false)
{
  fLEDModel = 0;
  fLETModel = 0;
  fLEAModel = 0;
  fTripathi = 0; 
  fTripathiLight = 0;
  fShen = 0;
  fIonH = 0;
  emax     = 20.*GeV;
  emaxLHEP = 1.*TeV;
  eminBIC  = 0.*MeV;
  SetPhysicsType(bIons);
  if(verbose > 1) G4cout << "### G4IonBinaryCascadePhysics" << G4endl;
}

G4IonBinaryCascadePhysics::~G4IonBinaryCascadePhysics()
{
  if(wasActivated) {
    delete fTripathi;
    delete fTripathiLight;
    delete fShen;
    delete fIonH;
    delete fLEDModel;
    delete fLETModel;
    delete fLEAModel;
    G4int i;
    G4int n = p_list.size();
    for(i=0; i<n; i++) {delete p_list[i];}
    n = model_list.size();
    for(i=0; i<n; i++) {delete model_list[i];}
  }
}

void G4IonBinaryCascadePhysics::ConstructProcess()
{
  if(wasActivated) { return; }
  wasActivated = true;

  G4BinaryLightIonReaction* fBC= new G4BinaryLightIonReaction();
  model_list.push_back(fBC);
  fShen = new G4IonsShenCrossSection;
  //fTripathi = new G4TripathiCrossSection;
  //fTripathiLight = new G4TripathiLightCrossSection;
  fIonH = new G4IonProtonCrossSection;

  fLEDModel = new G4LEDeuteronInelastic();
  fLETModel = new G4LETritonInelastic();
  fLEAModel = new G4LEAlphaInelastic();

  AddProcess("dInelastic", G4Deuteron::Deuteron(), fBC, fLEDModel);
  AddProcess("tInelastic",G4Triton::Triton(),  fBC, fLETModel);
  AddProcess("He3Inelastic",G4He3::He3(),  fBC, 0);
  AddProcess("alphaInelastic", G4Alpha::Alpha(),  fBC, fLEAModel);
  AddProcess("ionInelastic",G4GenericIon::GenericIon(),  fBC, 0);

}

void G4IonBinaryCascadePhysics::AddProcess(const G4String& name,
					   G4ParticleDefinition* p, 
					   G4HadronicInteraction* hmodel,
					   G4HadronicInteraction* lmodel)
{
  G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess(name, p);
  p_list.push_back(hadi);
  G4ProcessManager* pManager = p->GetProcessManager();
  pManager->AddDiscreteProcess(hadi);
  hadi->AddDataSet(fShen);
  //hadi->AddDataSet(fTripathi);
  //hadi->AddDataSet(fTripathiLight);
  if(p == G4GenericIon::GenericIon()) { hadi->AddDataSet(fIonH); }
  hmodel->SetMinEnergy(eminBIC);
  hmodel->SetMaxEnergy(emax);
  hadi->RegisterMe(hmodel);
  if(lmodel) {
    lmodel->SetMinEnergy(emax - MeV);
    lmodel->SetMaxEnergy(emaxLHEP);
    hadi->RegisterMe(lmodel);
  }  
  if(verbose > 1) {
    G4cout << "Register " << hadi->GetProcessName()
	   << " for " << p->GetParticleName()
	   << " Binary Cascade for E(MeV)= " << eminBIC << " - " << emax;
    if(lmodel) {
      G4cout  << " LHEP for E(MeV)= " << emax-MeV << " - " << emaxLHEP;
    }
    G4cout << G4endl;
  }
}

void G4IonBinaryCascadePhysics::ConstructParticle()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}
