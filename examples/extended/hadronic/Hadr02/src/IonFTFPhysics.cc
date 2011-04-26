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
// Class:    IonFTFPhysics
//
// Author:      A.Ivanchenko 02.03.2011
//
//
// Modified:
//
// ------------------------------------------------------------
// 

#include "IonFTFPhysics.hh"
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

#include "G4PreCompoundModel.hh"
#include "G4FTFBuilder.hh"
#include "G4HadronicInteraction.hh"
#include "G4PreCompoundModel.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IonFTFPhysics::IonFTFPhysics(G4bool val)
  : G4VHadronPhysics("ionInelasticFTF"),theIonBC(0),isBinary(val)
{
  fTripathi = fTripathiLight = fShen = fIonH = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IonFTFPhysics::~IonFTFPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IonFTFPhysics::ConstructProcess()
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

    fShen = new G4IonsShenCrossSection;
    fTripathi = new G4TripathiCrossSection;
    fTripathiLight = new G4TripathiLightCrossSection;
    fIonH = new G4IonProtonCrossSection;

    fShen->SetMaxKinEnergy(emax);
    fTripathi->SetMaxKinEnergy(emax);
    fTripathiLight->SetMaxKinEnergy(emax);
    fIonH->SetMaxKinEnergy(emax);    
  }
  G4PreCompoundModel* thePreCompound = 
    new G4PreCompoundModel(new G4ExcitationHandler());
  G4HadronicInteraction* theFTFP =
    BuildModel(new G4FTFBuilder("FTFP",thePreCompound),5.0*GeV,emax);
  //G4FTF2_5CrossSection *dpmXS = new G4FTF2_5CrossSection;
  G4cout << " IonFTFPhysics::ConstructProcess() FTFP+Binary for ions " 
	 << G4endl;

  /*
  // deuteron
  particle = G4Deuteron::Deuteron();
  G4cout << "IonFTFPhysics::ConstructProcess for " 
	 << particle->GetParticleName() << G4endl;
  pmanager = particle->GetProcessManager();
  G4HadronicProcess* hp = FindInelasticProcess(particle);
  hp->RegisterMe(theFTFP);
  pmanager->AddDiscreteProcess(hp);

  // triton
  particle = G4Triton::Triton();
  G4cout << "IonFTFPhysics::ConstructProcess for " 
	 << particle->GetParticleName() << G4endl;
  pmanager = particle->GetProcessManager();
  hp = FindInelasticProcess(particle);
  hp->RegisterMe(theFTFP);
  pmanager->AddDiscreteProcess(hp);

  // alpha
  particle = G4Alpha::Alpha();
  G4cout << "IonFTFPhysics::ConstructProcess for " 
	 << particle->GetParticleName() << G4endl;
  pmanager = particle->GetProcessManager();
  hp = FindInelasticProcess(particle);
  hp->RegisterMe(theFTFP);
  pmanager->AddDiscreteProcess(hp);
  */

  // He3
  particle = G4He3::He3();
  G4cout << "IonFTFPhysics::ConstructProcess for " 
	 << particle->GetParticleName() << G4endl;
  pmanager = particle->GetProcessManager();
  hp = FindInelasticProcess(particle);
  //hp->AddDataSet(dpmXS);
  if(!isBinary) { 
    hp->RegisterMe(theIonBC); 
    hp->AddDataSet(fShen);
    hp->AddDataSet(fTripathi);
    hp->AddDataSet(fTripathiLight);
    hp->AddDataSet(fIonH); 
  }
  hp->RegisterMe(theFTFP);
  pmanager->AddDiscreteProcess(hp);

  // GenericIon
  particle = G4GenericIon::GenericIon();
  G4cout << "IonFTFPhysics::ConstructProcess for " 
	 << particle->GetParticleName() << G4endl;
  pmanager = particle->GetProcessManager();
  hp = FindInelasticProcess(particle);
  G4cout << "IonFTFPhysics: FTF is registered for GenericIon" << G4endl;
  //hp->AddDataSet(dpmXS);
  if(!isBinary) { 
    hp->RegisterMe(theIonBC); 
    hp->AddDataSet(fShen);
    hp->AddDataSet(fTripathi);
    hp->AddDataSet(fTripathiLight);
    hp->AddDataSet(fIonH); 
  }
  hp->RegisterMe(theFTFP);
  pmanager->AddDiscreteProcess(hp);

  G4cout << "IonFTFPhysics::ConstructProcess done! " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



