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
// $Id: PhysListIonBinaryCascade.cc,v 1.5 2005/03/16 13:10:37 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysListIonBinaryCascade.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListIonBinaryCascade::PhysListIonBinaryCascade(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListIonBinaryCascade::~PhysListIonBinaryCascade()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListIonBinaryCascade::ConstructProcess()
{
  // Binary Cascade
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* pmanager = 0;
  G4BinaryLightIonReaction* theBC = new G4BinaryLightIonReaction();
  theBC->SetMinEnergy(80*MeV);
  theBC->SetMaxEnergy(20*GeV);

  G4TripathiCrossSection * TripathiCrossSection= new G4TripathiCrossSection;
  G4IonsShenCrossSection * aShen = new G4IonsShenCrossSection;

  // deuteron
  particle = G4Deuteron::Deuteron();
  pmanager = particle->GetProcessManager();
  G4LEDeuteronInelastic* theDIModel = new G4LEDeuteronInelastic;
  theDIModel->SetMaxEnergy(100*MeV);
  theIPdeuteron.AddDataSet(TripathiCrossSection);
  theIPdeuteron.AddDataSet(aShen);
  theIPdeuteron.RegisterMe(theDIModel);
  theIPdeuteron.RegisterMe(theBC);
  pmanager->AddDiscreteProcess(&theIPdeuteron);

  // triton
  particle = G4Triton::Triton();
  pmanager = particle->GetProcessManager();
  G4LETritonInelastic* theTIModel = new G4LETritonInelastic;
  theTIModel->SetMaxEnergy(100*MeV);
  theIPtriton.AddDataSet(TripathiCrossSection);
  theIPtriton.AddDataSet(aShen);
  theIPtriton.RegisterMe(theTIModel);
  theIPtriton.RegisterMe(theBC);
  pmanager->AddDiscreteProcess(&theIPtriton);

  // alpha
  particle = G4Alpha::Alpha();
  pmanager = particle->GetProcessManager();
  G4LEAlphaInelastic* theAIModel = new G4LEAlphaInelastic;
  theAIModel->SetMaxEnergy(100*MeV);
  theIPalpha.AddDataSet(TripathiCrossSection);
  theIPalpha.AddDataSet(aShen);
  theIPalpha.RegisterMe(theAIModel);
  theIPalpha.RegisterMe(theBC);
  pmanager->AddDiscreteProcess(&theIPalpha);

  // GenericIon
  particle = G4GenericIon::GenericIon();
  pmanager = particle->GetProcessManager();
  G4HadronInelasticProcess* theIPGenericIon =
      new G4HadronInelasticProcess("IonInelastic",particle);
  theIPGenericIon->AddDataSet(TripathiCrossSection);
  theIPGenericIon->AddDataSet(aShen);
  G4BinaryLightIonReaction * theGenIonBC= new G4BinaryLightIonReaction;
  theGenIonBC->SetMinEnergy(0*MeV);
  theGenIonBC->SetMaxEnergy(10*GeV);
  theIPGenericIon->RegisterMe(theGenIonBC);
  pmanager->AddDiscreteProcess(theIPGenericIon);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

