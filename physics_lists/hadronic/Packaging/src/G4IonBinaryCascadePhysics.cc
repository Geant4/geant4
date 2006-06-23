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
// $Id: G4IonBinaryCascadePhysics.cc,v 1.3 2006-06-23 08:23:03 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4IonBinaryCascadePhysics
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
// 23.06.06 V.Ivanchenko set emaxLHEP=1 TeV
//
//----------------------------------------------------------------------------
//

#include "G4IonBinaryCascadePhysics.hh"

#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"

#include "G4BinaryLightIonReaction.hh"
#include "G4TripathiCrossSection.hh"
#include "G4TripathiLightCrossSection.hh"
#include "G4IonsShenCrossSection.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

// Nuclei
#include "G4IonConstructor.hh"

G4IonBinaryCascadePhysics::G4IonBinaryCascadePhysics(const G4String& name, 
						     G4int verb)
  :  G4VPhysicsConstructor(name), verbose(verb), wasActivated(false)
{
  emax     = 20.*GeV;
  emaxLHEP = 1.*TeV;
  eminBIC  = 0.*MeV;
  if(verbose > 1) G4cout << "### G4IonBinaryCascadePhysics" << G4endl;
}

G4IonBinaryCascadePhysics::~G4IonBinaryCascadePhysics()
{
  if(wasActivated) {
    delete fTripathi;
    delete fTripathiLight;
    delete fShen;
    G4int i;
    G4int n = p_list.size();
    for(i=0; i<n; i++) {delete p_list[i];}
    n = model_list.size();
    for(i=0; i<n; i++) {delete model_list[i];}
  }
}

void G4IonBinaryCascadePhysics::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  G4BinaryLightIonReaction* fBC= new G4BinaryLightIonReaction();
  model_list.push_back(fBC);
  fShen = new G4IonsShenCrossSection;
  fTripathi = new G4TripathiCrossSection;
  fTripathiLight = new G4TripathiLightCrossSection;

  //    new G4LEDeuteronInelastic;

  AddProcess("dInelastic", G4Deuteron::Deuteron(), fBC, 0);
  AddProcess("tInelastic",G4Triton::Triton(),  fBC, 0);
  AddProcess("He3Inelastic",G4He3::He3(),  fBC, 0);
  AddProcess("alphaInelastic", G4Alpha::Alpha(),  fBC, 0);
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
  hadi->AddDataSet(fTripathi);
  hadi->AddDataSet(fTripathiLight);
  hadi->AddDataSet(fShen);
  hadi->RegisterMe(hmodel);
  hmodel->SetMinEnergy(eminBIC);
  hmodel->SetMaxEnergy(emax);
  if(lmodel) {
    hadi->RegisterMe(lmodel);
    hmodel->SetMinEnergy(emax - MeV);
    hmodel->SetMaxEnergy(emaxLHEP);
  }  
  if(verbose > 1) {
    G4cout << "Register " << hadi->GetProcessName()
	   << " for " << p->GetParticleName()
	   << " Binary Cascade for E(MeV)= " << eminBIC << " - " << emax 
	   <<G4endl;
  }
}

void G4IonBinaryCascadePhysics::ConstructParticle()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}
