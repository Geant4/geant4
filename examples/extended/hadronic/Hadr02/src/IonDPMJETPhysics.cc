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

#include "G4BuilderType.hh"

#ifdef G4_USE_DPMJET
#include "G4DPMJET2_5Model.hh"
#include "G4DPMJET2_5Interface.hh"
#include "G4DPMJET2_5CrossSection.hh"
#endif

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IonDPMJETPhysics::IonDPMJETPhysics(G4bool val)
  : G4VHadronPhysics("ionInelasticDPMJET"),theIonBC(0),theDPM(0),
    useDPMJETXS(val)
{
  fTripathi = fTripathiLight = fShen = fIonH = 0;
  SetPhysicsType(bIons);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IonDPMJETPhysics::~IonDPMJETPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IonDPMJETPhysics::ConstructProcess()
{
  G4double emax = 1000.*TeV;

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

#ifdef G4_USE_DPMJET
  theDPM = new G4DPMJET2_5Model();
  theDPM->SetMinEnergy(5*GeV);
  theDPM->SetMaxEnergy(emax);
  //G4DPMJET2_5CrossSection *dpmXS = new G4DPMJET2_5CrossSection;
#endif

  AddProcess("dInelastic", G4Deuteron::Deuteron(),false);
  AddProcess("tInelastic",G4Triton::Triton(),false);
  AddProcess("He3Inelastic",G4He3::He3(),true);
  AddProcess("alphaInelastic", G4Alpha::Alpha(),true);
  AddProcess("ionInelastic",G4GenericIon::GenericIon(),true);

  G4cout << "IonDPMJETPhysics::ConstructProcess done! " << G4endl;
}

void IonDPMJETPhysics::AddProcess(const G4String& name,
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
#ifdef G4_USE_DPMJET
  hadi->RegisterMe(theDPM); 
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



