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
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"

#include "G4BinaryLightIonReaction.hh"
#include "G4QMDReaction.hh"

//#include "G4TripathiCrossSection.hh"
//#include "G4TripathiLightCrossSection.hh"
//#include "G4IonsShenCrossSection.hh"

#include "G4GGNuclNuclCrossSection.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

// Nuclei
#include "G4IonConstructor.hh"
#include "G4BuilderType.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4IonQMDPhysics);

G4IonQMDPhysics::G4IonQMDPhysics(G4int ver)
  :  G4VPhysicsConstructor("IonQMD"), verbose(ver), wasActivated(false)
{
  fLEDModel = 0;
  fLETModel = 0;
  fLEAModel = 0;
//  fTripathi = 0;
//  fTripathiLight = 0;
//  fShen = 0;
  eminBIC  = 0.*MeV;
  eminQMD  = 100.*MeV;
  emaxQMD  = 10.*GeV;
  emaxLHEP = 1.*TeV;
  overlap  = 10*MeV;
  SetPhysicsType(bIons);
  if(verbose > 1) G4cout << "### G4IonQMDPhysics" << G4endl;
}

G4IonQMDPhysics::G4IonQMDPhysics(const G4String& name, 
						     G4int ver)
  :  G4VPhysicsConstructor(name), verbose(ver), wasActivated(false)
{
  fLEDModel = 0;
  fLETModel = 0;
  fLEAModel = 0;
//  fTripathi = 0;
//  fTripathiLight = 0;
//  fShen = 0;
  eminBIC  = 0.*MeV;
  eminQMD  = 100.*MeV;
  emaxQMD  = 10.*GeV;
  emaxLHEP = 1.*TeV;
  overlap  = 10*MeV;
  SetPhysicsType(bIons);
  if(verbose > 1) G4cout << "### G4IonQMDPhysics" << G4endl;
}

G4IonQMDPhysics::~G4IonQMDPhysics()
{
  if(wasActivated) {
//    delete fTripathi;
//    delete fTripathiLight;
//    delete fShen;
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

void G4IonQMDPhysics::ConstructProcess()
{
  if(wasActivated) return;
  wasActivated = true;

  G4BinaryLightIonReaction* fBC= new G4BinaryLightIonReaction();
  model_list.push_back(fBC);
  G4QMDReaction* fQMD= new G4QMDReaction();
  model_list.push_back(fQMD);
  
//  fShen = new G4IonsShenCrossSection;
//  fTripathi = new G4TripathiCrossSection
//    fTripathiLight = new G4TripathiLightCrossSection;

  fLEDModel = new G4LEDeuteronInelastic();
  fLETModel = new G4LETritonInelastic();
  fLEAModel = new G4LEAlphaInelastic();

  AddProcess("dInelastic", G4Deuteron::Deuteron(), fBC, fQMD, fLEDModel );
  AddProcess("tInelastic",G4Triton::Triton(),      fBC, fQMD, fLETModel );
  AddProcess("He3Inelastic",G4He3::He3(),          fBC, fQMD, 0 );
  AddProcess("alphaInelastic", G4Alpha::Alpha(),   fBC, fQMD, fLEAModel );
  AddProcess("ionInelastic",G4GenericIon::GenericIon(), fBC, fQMD, 0);

}

void G4IonQMDPhysics::AddProcess(const G4String& name,
					   G4ParticleDefinition* p, 
					   G4BinaryLightIonReaction* BIC,
					   G4QMDReaction* QMD,
					   G4HadronicInteraction* LHEP)
{
  G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess(name, p);
  p_list.push_back(hadi);
  G4ProcessManager* pManager = p->GetProcessManager();
  pManager->AddDiscreteProcess(hadi);
  
//  hadi->AddDataSet(fShen);
//  hadi->AddDataSet(fTripathi);
//  hadi->AddDataSet(fTripathiLight);

  hadi->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4GGNuclNuclCrossSection::Default_Name()));
    
  BIC->SetMinEnergy(eminBIC);
  BIC->SetMaxEnergy(emaxQMD);  //reset when QMD is present 
  hadi->RegisterMe(BIC);

  if(QMD) {
    QMD->SetMinEnergy(eminQMD);
    BIC->SetMaxEnergy(eminQMD+overlap);
    QMD->SetMaxEnergy(emaxQMD);
    hadi->RegisterMe(QMD);
  }  

  if(LHEP) {
    LHEP->SetMinEnergy(emaxQMD - overlap);
    LHEP->SetMaxEnergy(emaxLHEP);
    hadi->RegisterMe(LHEP);
  }  
  if(verbose > 1) {
    G4cout << "Register " << hadi->GetProcessName()
	   << " for " << p->GetParticleName() << G4endl
	   << "       Binary Cascade for E(MeV)= " << eminBIC << " - " 
	   << (QMD==0 ? emaxQMD : (eminQMD-overlap)) ;
    if(QMD) {
      G4cout  << G4endl <<"       QMD for E(MeV)= " << eminQMD << " - " << emaxQMD;
    }
    if(LHEP) {
      G4cout << G4endl<< "       LHEP for E(MeV)= " << emaxQMD-overlap << " - " << emaxLHEP;
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
