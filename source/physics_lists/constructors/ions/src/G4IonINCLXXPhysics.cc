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
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"

#include "G4INCLXXInterface.hh"
#include "G4TripathiCrossSection.hh"
#include "G4TripathiLightCrossSection.hh"
#include "G4IonsShenCrossSection.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

// Nuclei
#include "G4IonConstructor.hh"
#include "G4BuilderType.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4IonINCLXXPhysics);


G4IonINCLXXPhysics::G4IonINCLXXPhysics(G4int ver) :
  G4VPhysicsConstructor("IonINCLXX"),
  fINCLXXIons(NULL),
  fTripathi(NULL),
  fTripathiLight(NULL),
  fShen(NULL),
  fLEDModel(NULL),
  fLETModel(NULL),
  fLEAModel(NULL),
  verbose(ver), wasActivated(false)
{
  // INCLXX light ion maximum energy is 3.0 GeV/nucleon
  emax_d     = 2 * 3.0 * GeV;
  emax_t     = 3 * 3.0 * GeV;
  emax_he3   = 3 * 3.0 * GeV;
  emax_alpha = 4 * 3.0 * GeV;
  emax       = 16 * 3.0 * GeV;
  emaxLHEP   = 1.*TeV;
  emin       = 0.*MeV;
  SetPhysicsType(bIons);
  if(verbose > 1) G4cout << "### G4IonINCLXXPhysics" << G4endl;
}

G4IonINCLXXPhysics::G4IonINCLXXPhysics(const G4String& name, 
						     G4int ver)
  :  G4VPhysicsConstructor(name),
  fINCLXXIons(NULL),
  fTripathi(NULL),
  fTripathiLight(NULL),
  fShen(NULL),
  fLEDModel(NULL),
  fLETModel(NULL),
  fLEAModel(NULL),
  verbose(ver), wasActivated(false)
{
  // INCLXX light ion maximum energy is 3.0 GeV/nucleon
  emax_d     = 2 * 3.0 * GeV;
  emax_t     = 3 * 3.0 * GeV;
  emax_he3   = 3 * 3.0 * GeV;
  emax_alpha = 4 * 3.0 * GeV;
  emax       = 16 * 3.0 * GeV;
  emaxLHEP   = 1.*TeV;
  emin       = 0.*MeV;
  SetPhysicsType(bIons);
  if(verbose > 1) G4cout << "### G4IonINCLXXPhysics" << G4endl;
}

G4IonINCLXXPhysics::~G4IonINCLXXPhysics()
{
  if(wasActivated) {
    delete fTripathi;
    delete fTripathiLight;
    delete fShen;
    delete fLEDModel;
    delete fLETModel;
    delete fLEAModel;
    G4int i;
    G4int n = p_list.size();
    for(i=0; i<n; i++) {delete p_list[i];}
    delete fINCLXXIons;
  }
}

void G4IonINCLXXPhysics::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  fINCLXXIons= new G4INCLXXInterface();
  fShen = new G4IonsShenCrossSection;
  fTripathi = new G4TripathiCrossSection;
  fTripathiLight = new G4TripathiLightCrossSection;

  fLEDModel = new G4LEDeuteronInelastic();
  fLETModel = new G4LETritonInelastic();
  fLEAModel = new G4LEAlphaInelastic();

  AddProcess("dInelastic", G4Deuteron::Deuteron(), fINCLXXIons, fLEDModel, emax_d);
  AddProcess("tInelastic",G4Triton::Triton(),  fINCLXXIons, fLETModel, emax_t);
  AddProcess("He3Inelastic",G4He3::He3(),  fINCLXXIons, 0, emax_he3);
  AddProcess("alphaInelastic", G4Alpha::Alpha(),  fINCLXXIons, fLEAModel, emax_alpha);
  AddProcess("ionInelastic",G4GenericIon::GenericIon(),  fINCLXXIons, 0, emax);
}

void G4IonINCLXXPhysics::AddProcess(const G4String& name,
					   G4ParticleDefinition* p, 
					   G4HadronicInteraction* hmodel,
					   G4HadronicInteraction* lmodel,
				      const G4double inclxxEnergyUpperLimit = 3.0 * GeV)
{
  G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess(name, p);
  p_list.push_back(hadi);
  G4ProcessManager* pManager = p->GetProcessManager();
  pManager->AddDiscreteProcess(hadi);
  hadi->AddDataSet(fShen);
  hadi->AddDataSet(fTripathi);
  hadi->AddDataSet(fTripathiLight);
  hmodel->SetMinEnergy(emin);
  hmodel->SetMaxEnergy(inclxxEnergyUpperLimit);
  hadi->RegisterMe(hmodel);
  if(lmodel) {
    lmodel->SetMinEnergy(inclxxEnergyUpperLimit - MeV);
    lmodel->SetMaxEnergy(emaxLHEP);
    hadi->RegisterMe(lmodel);
  }
  if(verbose > 1) {
    G4cout << "Register " << hadi->GetProcessName()
	   << " for " << p->GetParticleName()
	   << " INCLXX/G4DeexcitationHandler for E(MeV)= " << emin << " - " << inclxxEnergyUpperLimit;
    if(lmodel) {
      G4cout  << " LHEP for E(MeV)= " << inclxxEnergyUpperLimit-MeV << " - " << emaxLHEP;
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
