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
// $Id: PhysListIonBinaryCascade.cc,v 1.1 2003-11-19 10:16:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

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

#include "G4HadronInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"

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
  G4BinaryLightIonReaction* theBC = 0;

  // deuteron
  particle = G4Deuteron::Deuteron();
  pmanager = particle->GetProcessManager();
  theBC = new G4BinaryLightIonReaction();
  theIPdeuteron.RegisterMe(theBC);
  pmanager->AddDiscreteProcess(&theIPdeuteron);

  // triton
  particle = G4Triton::Triton();
  pmanager = particle->GetProcessManager();
  theBC = new G4BinaryLightIonReaction();
  theIPtriton.RegisterMe(theBC);
  pmanager->AddDiscreteProcess(&theIPtriton);

  // alpha
  particle = G4Alpha::Alpha();
  pmanager = particle->GetProcessManager();
  theBC = new G4BinaryLightIonReaction();
  theIPalpha.RegisterMe(theBC);
  pmanager->AddDiscreteProcess(&theIPalpha);

  // He3
  particle = G4He3::He3();
  pmanager = particle->GetProcessManager();
  G4HadronInelasticProcess* theIPHe3 =
      new G4HadronInelasticProcess("He3Inelastic",particle);
  theBC = new G4BinaryLightIonReaction();
  theIPHe3->RegisterMe(theBC);
  pmanager->AddDiscreteProcess(theIPHe3);

  // C12
  particle = IonC12::Ion();
  pmanager = particle->GetProcessManager();
  G4HadronInelasticProcess* theIPIonC12 =
      new G4HadronInelasticProcess("IonC12Inelastic",particle);
  theBC = new G4BinaryLightIonReaction();
  theIPIonC12->RegisterMe(theBC);
  pmanager->AddDiscreteProcess(theIPIonC12);

  // GenericIon
  particle = G4GenericIon::GenericIon();
  pmanager = particle->GetProcessManager();
  G4HadronInelasticProcess* theIPGenericIon =
      new G4HadronInelasticProcess("IonInelastic",particle);
  theBC = new G4BinaryLightIonReaction();
  theIPGenericIon->RegisterMe(theBC);
  pmanager->AddDiscreteProcess(theIPGenericIon);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

