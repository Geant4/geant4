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
// $Id: G4IonInclAblaPhysics.cc,v 1.2 2010-06-03 15:03:53 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4IonInclAblaPhysics
//
// Author:      P. Kaitaniemi
//
// Modified:
// 23.06.06 V.Ivanchenko set emaxLHEP=1 TeV
// 24.06.06 V.Ivanchenko fix typo
//
//----------------------------------------------------------------------------
//

#include "G4IonInclAblaPhysics.hh"

#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"

#include "G4InclAblaLightIonInterface.hh"
#include "G4TripathiCrossSection.hh"
#include "G4TripathiLightCrossSection.hh"
#include "G4IonsShenCrossSection.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

// Nuclei
#include "G4IonConstructor.hh"

G4IonInclAblaPhysics::G4IonInclAblaPhysics(G4int ver)
  :  G4VPhysicsConstructor("IonInclAbla"), verbose(ver), wasActivated(false)
{
  // INCL/ABLA light ion maximum energy is 3.0 GeV/nucleon
  emax_d     = 2 * 3.0 * GeV;
  emax_t     = 3 * 3.0 * GeV;
  emax_he3   = 3 * 3.0 * GeV;
  emax_alpha = 4 * 3.0 * GeV;
  emaxLHEP   = 1.*TeV;
  emin       = 0.*MeV;
  if(verbose > 1) G4cout << "### G4IonInclAblaPhysics" << G4endl;
}

G4IonInclAblaPhysics::G4IonInclAblaPhysics(const G4String& name, 
						     G4int ver)
  :  G4VPhysicsConstructor(name), verbose(ver), wasActivated(false)
{
  // INCL/ABLA light ion maximum energy is 3.0 GeV/nucleon
  emax_d     = 2 * 3.0 * GeV;
  emax_t     = 3 * 3.0 * GeV;
  emax_he3   = 3 * 3.0 * GeV;
  emax_alpha = 4 * 3.0 * GeV;
  emaxLHEP   = 1.*TeV;
  emin       = 0.*MeV;
  if(verbose > 1) G4cout << "### G4IonInclAblaPhysics" << G4endl;
}

G4IonInclAblaPhysics::~G4IonInclAblaPhysics()
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
    n = model_list.size();
    for(i=0; i<n; i++) {delete model_list[i];}
  }
}

void G4IonInclAblaPhysics::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  G4InclAblaLightIonInterface* fInclAblaIons= new G4InclAblaLightIonInterface();
  model_list.push_back(fInclAblaIons);
  fShen = new G4IonsShenCrossSection;
  fTripathi = new G4TripathiCrossSection;
  fTripathiLight = new G4TripathiLightCrossSection;

  fLEDModel = new G4LEDeuteronInelastic();
  fLETModel = new G4LETritonInelastic();
  fLEAModel = new G4LEAlphaInelastic();

  AddProcess("dInelastic", G4Deuteron::Deuteron(), fInclAblaIons, fLEDModel, emax_d);
  AddProcess("tInelastic",G4Triton::Triton(),  fInclAblaIons, fLETModel, emax_t);
  AddProcess("He3Inelastic",G4He3::He3(),  fInclAblaIons, 0, emax_he3);
  AddProcess("alphaInelastic", G4Alpha::Alpha(),  fInclAblaIons, fLEAModel, emax_alpha);
  // Support for light ions heavier than Alpha will be included in a future release of INCL/ABLA
  //  AddProcess("ionInelastic",G4GenericIon::GenericIon(),  fInclAblaIons, 0);
}

void G4IonInclAblaPhysics::AddProcess(const G4String& name,
					   G4ParticleDefinition* p, 
					   G4HadronicInteraction* hmodel,
					   G4HadronicInteraction* lmodel,
				      const G4double inclEnergyUpperLimit = 3.0 * GeV)
{
  G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess(name, p);
  p_list.push_back(hadi);
  G4ProcessManager* pManager = p->GetProcessManager();
  pManager->AddDiscreteProcess(hadi);
  hadi->AddDataSet(fShen);
  hadi->AddDataSet(fTripathi);
  hadi->AddDataSet(fTripathiLight);
  hmodel->SetMinEnergy(emin);
  hmodel->SetMaxEnergy(inclEnergyUpperLimit);
  hadi->RegisterMe(hmodel);
  if(lmodel) {
    lmodel->SetMinEnergy(inclEnergyUpperLimit - MeV);
    lmodel->SetMaxEnergy(emaxLHEP);
    hadi->RegisterMe(lmodel);
  }  
  if(verbose > 1) {
    G4cout << "Register " << hadi->GetProcessName()
	   << " for " << p->GetParticleName()
	   << " INCL/ABLA for E(MeV)= " << emin << " - " << inclEnergyUpperLimit;
    if(lmodel) {
      G4cout  << " LHEP for E(MeV)= " << inclEnergyUpperLimit-MeV << " - " << emaxLHEP;
    }
    G4cout << G4endl;
  }
}

void G4IonInclAblaPhysics::ConstructParticle()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}
