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
// $Id: IonABPhysics.cc,v 1.1 2006/10/28 16:00:25 vnivanch Exp $
// GEANT4 tag $Name: gras-02-05-02 $
//
//---------------------------------------------------------------------------
//
// Class:    IonDPMJETPhysics
//
// Author:      A.Ivanchenko 26.08.2010
//
// This class was designed under ESA contracts
// 
// Customer:     
// Contract:            
//
//
// Modified:
//
// ------------------------------------------------------------
// 

#include "IonDPMJETPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4TripathiCrossSection.hh"
#include "G4TripathiLightCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4IonProtonCrossSection.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"

#include "G4PreCompoundModel.hh"

#ifdef G4_USE_DPMJET
#include "G4DPMJET2_5Model.hh"
#include "G4DPMJET2_5Interface.hh"
#include "G4DPMJET2_5CrossSection.hh"
#endif

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IonDPMJETPhysics::IonDPMJETPhysics(G4bool val)
  : G4VHadronPhysics("ionInelasticDPMJET"),theIonBC(0),theDPM(0),isBinary(val)
{
  fTripathi = fTripathiLight = fShen = fIonH = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IonDPMJETPhysics::~IonDPMJETPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IonDPMJETPhysics::ConstructProcess()
{
  // Binary Cascade
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* pmanager = 0;
  G4HadronicProcess* hp = 0;

  G4double emax = 1000.*TeV;

  if(!isBinary) {
    theIonBC = new G4BinaryLightIonReaction();
    theIonBC->SetMinEnergy(0.0);
    theIonBC->SetMaxEnergy(6*GeV);

    fShen = new G4IonsShenCrossSection();
    fTripathi = new G4TripathiCrossSection();
    fTripathiLight = new G4TripathiLightCrossSection();
    fIonH = new G4IonProtonCrossSection();

    fShen->SetMaxKinEnergy(emax);
    fTripathi->SetMaxKinEnergy(emax);
    fTripathiLight->SetMaxKinEnergy(emax);
    fIonH->SetMaxKinEnergy(emax);    
  }

#ifdef G4_USE_DPMJET
  theDPM = new G4DPMJET2_5Model();
  theDPM->SetMinEnergy(5*GeV);
  theDPM->SetMaxEnergy(emax);
  //G4DPMJET2_5CrossSection *dpmXS = new G4DPMJET2_5CrossSection;
  G4cout << " IonDPMJETPhysics::ConstructProcess() DPMJET-> " 
	 << theDPM << G4endl;
#endif


  // not possible to use for d,t,alpha
  /*
  // deuteron
  particle = G4Deuteron::Deuteron();
  G4cout << "IonDPMJETPhysics::ConstructProcess for " 
	 << particle->GetParticleName() << G4endl;
  pmanager = particle->GetProcessManager();
  hp = FindInelasticProcess(particle);
#ifdef G4_USE_DPMJET
  //hp->AddDataSet(dpmXS);
  hp->RegisterMe(theDPM);
#endif
  pmanager->AddDiscreteProcess(hp);

  // triton
  particle = G4Triton::Triton();
  G4cout << "IonDPMJETPhysics::ConstructProcess for " 
	 << particle->GetParticleName() << G4endl;
  pmanager = particle->GetProcessManager();
  hp = FindInelasticProcess(particle);
#ifdef G4_USE_DPMJET
  //hp->AddDataSet(dpmXS);
  hp->RegisterMe(theDPM);
#endif
  pmanager->AddDiscreteProcess(hp);

  // alpha
  particle = G4Alpha::Alpha();
  G4cout << "IonDPMJETPhysics::ConstructProcess for " 
	 << particle->GetParticleName() << G4endl;
  pmanager = particle->GetProcessManager();
  hp = FindInelasticProcess(particle);
#ifdef G4_USE_DPMJET
  //hp->AddDataSet(dpmXS);
  hp->RegisterMe(theDPM);
#endif
  pmanager->AddDiscreteProcess(hp);
  */

  // He3
  particle = G4He3::He3();
  G4cout << "IonDPMJETPhysics::ConstructProcess for " 
	 << particle->GetParticleName() << G4endl;
  pmanager = particle->GetProcessManager();
  hp = FindInelasticProcess(particle);
  if(!isBinary) { 
    hp->RegisterMe(theIonBC); 
    hp->AddDataSet(fShen);
    hp->AddDataSet(fTripathi);
    hp->AddDataSet(fTripathiLight);
    hp->AddDataSet(fIonH); 
  }
#ifdef G4_USE_DPMJET
  //hp->AddDataSet(dpmXS);
  hp->RegisterMe(theDPM);
#endif
  pmanager->AddDiscreteProcess(hp);

  // GenericIon
  particle = G4GenericIon::GenericIon();
  G4cout << "IonDPMJETPhysics::ConstructProcess for " 
	 << particle->GetParticleName() << G4endl;
  pmanager = particle->GetProcessManager();
  hp = FindInelasticProcess(particle);
  if(!isBinary) { 
    hp->RegisterMe(theIonBC); 
    hp->AddDataSet(fShen);
    hp->AddDataSet(fTripathi);
    hp->AddDataSet(fTripathiLight);
    hp->AddDataSet(fIonH); 
  }
#ifdef G4_USE_DPMJET
  //hp->AddDataSet(dpmXS);
  hp->RegisterMe(theDPM);
  G4cout << "IonDPMJETPhysics: DPMJET is registered for GenericIon" << G4endl;
#endif
  pmanager->AddDiscreteProcess(hp);

  G4cout << "IonDPMJETPhysics::ConstructProcess done! " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



