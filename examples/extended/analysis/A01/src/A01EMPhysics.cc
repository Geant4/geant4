// $Id: A01EMPhysics.cc,v 1.1 2002-11-13 07:23:01 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//


#include "A01EMPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/iomanip"


A01EMPhysics::A01EMPhysics(const G4String& name)
               :  G4VPhysicsConstructor(name)
{
}

A01EMPhysics::~A01EMPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4Gamma.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"

void A01EMPhysics::ConstructParticle()
{
  // gamma
  G4Gamma::GammaDefinition();

  // electron
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
}


#include "G4ProcessManager.hh"


void A01EMPhysics::ConstructProcess()
{
  G4ProcessManager * pManager = 0;

  // Gamma Physics
  pManager = G4Gamma::Gamma()->GetProcessManager();
  pManager->AddDiscreteProcess(&thePhotoEffect);
  pManager->AddDiscreteProcess(&theComptonEffect);
  pManager->AddDiscreteProcess(&thePairProduction);

  // Electron Physics
  pManager = G4Electron::Electron()->GetProcessManager();
   // add processes
  pManager->AddDiscreteProcess(&theElectronBremsStrahlung);

  pManager->AddProcess(&theElectronIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theElectronMultipleScattering);
  pManager->SetProcessOrdering(&theElectronMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&theElectronMultipleScattering, idxPostStep,  1);

  //Positron Physics
  pManager = G4Positron::Positron()->GetProcessManager();
  // add processes
  pManager->AddDiscreteProcess(&thePositronBremsStrahlung);

  pManager->AddDiscreteProcess(&theAnnihilation);

  pManager->AddRestProcess(&theAnnihilation);

  pManager->AddProcess(&thePositronIonisation, ordInActive,2, 2);

  pManager->AddProcess(&thePositronMultipleScattering);
  pManager->SetProcessOrdering(&thePositronMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&thePositronMultipleScattering, idxPostStep,  1);

}



