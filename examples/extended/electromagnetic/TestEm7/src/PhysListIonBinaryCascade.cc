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
// $Id: PhysListIonBinaryCascade.cc,v 1.2 2003/12/05 11:17:16 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-00 $

#include "PhysListIonBinaryCascade.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "IonC12.hh"
#include "G4GenericIon.hh"

#include "G4ExcitationHandler.hh"
#include "G4Evaporation.hh"
#include "G4FermiBreakUp.hh"
#include "G4StatMF.hh"
#include "G4GeneratorPrecompoundInterface.hh"

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
  /*
  // all models for treatment of thermal nucleus
  G4Evaporation * theEvaporation = new G4Evaporation;
  G4FermiBreakUp * theFermiBreakUp = new G4FermiBreakUp;
  G4StatMF * theMF = new G4StatMF;

  // Evaporation logic
  G4ExcitationHandler * theHandler = new G4ExcitationHandler;
  theHandler->SetEvaporation(theEvaporation);
  theHandler->SetFermiModel(theFermiBreakUp);
  theHandler->SetMultiFragmentation(theMF);
  theHandler->SetMaxAandZForFermiBreakUp(12, 6);
  theHandler->SetMinEForMultiFrag(3*MeV);

  // Pre equilibrium stage
  G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel(theHandler);

  // a no-cascade generator-precompound interaface
  G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface;
             theCascade->SetDeExcitation(thePreEquilib);
  */
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

  // C12
  particle = IonC12::Ion();
  pmanager = particle->GetProcessManager();
  G4HadronInelasticProcess* theIPIonC12 =
      new G4HadronInelasticProcess("IonC12Inelastic",particle);
  theIPIonC12->AddDataSet(TripathiCrossSection);
  theIPIonC12->AddDataSet(aShen);
  theIPIonC12->RegisterMe(theGenIonBC);
  pmanager->AddDiscreteProcess(theIPIonC12);

  // He3
  particle = G4He3::He3();
  pmanager = particle->GetProcessManager();
  G4HadronInelasticProcess* theIPHe3 =
      new G4HadronInelasticProcess("He3Inelastic",particle);
  theIPHe3->AddDataSet(TripathiCrossSection);
  theIPHe3->AddDataSet(aShen);
  theIPHe3->RegisterMe(theGenIonBC);
  pmanager->AddDiscreteProcess(theIPHe3);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

