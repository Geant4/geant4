#include "G4EMBuilder.hh"

#include "globals.hh"
#include "G4ios.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"

G4EMBuilder::
G4EMBuilder(): wasActivated(false), synchOn(false), 
               gammNucOn(false)
{
  theT = new G4EMTailorer(this);
}

void G4EMBuilder::
Synch (G4String & newState)
{
  if(newState == "on")
  {
    if(synchOn)
    {
    }
    else if(wasActivated)
    {
      G4ProcessManager * pManager = 0;
      pManager = G4Electron::ElectronDefinition()->GetProcessManager();
      pManager->AddDiscreteProcess(&theElectronSynch);
      pManager = G4Positron::PositronDefinition()->GetProcessManager();
      pManager->AddDiscreteProcess(&thePositronSynch);
    }
    synchOn = true;
  }
  else if(synchOn)
  {
    synchOn = false;
    G4ProcessManager * pManager = 0;
    pManager = G4Electron::ElectronDefinition()->GetProcessManager();
    pManager->RemoveProcess(&theElectronSynch);
    pManager = G4Positron::PositronDefinition()->GetProcessManager();
    pManager->RemoveProcess(&thePositronSynch);
  }
}

void G4EMBuilder::
GammaNuclear (G4String & newState)
{
  if(newState == "on")
  {
    if(gammNucOn)
    {
    }
    else if(wasActivated)  
    {
      theGNPhysics.Build();
    }
    gammNucOn = true;
  }
  else if(gammNucOn)
  {
    G4Exception("Switching off of gamma nuclear currently not supported.");
  }
}

G4EMBuilder::
~G4EMBuilder() 
{
  if(wasActivated)
  {
  G4ProcessManager * pManager = 0;
  pManager = G4Gamma::Gamma()->GetProcessManager();
  if(pManager)
  {
  pManager->RemoveProcess(&thePhotoEffect);
  pManager->RemoveProcess(&theComptonEffect);
  pManager->RemoveProcess(&thePairProduction);
  }

  pManager = G4Electron::Electron()->GetProcessManager();
  if(pManager)
  {
  pManager->RemoveProcess(&theElectronBremsStrahlung);  
  pManager->RemoveProcess(&theElectronIonisation);
  pManager->RemoveProcess(&theElectronMultipleScattering);
  if(synchOn) pManager->RemoveProcess(&theElectronSynch);
  }

  pManager = G4Positron::Positron()->GetProcessManager();
  if(pManager)
  {
  pManager->RemoveProcess(&thePositronBremsStrahlung);
  pManager->RemoveProcess(&theAnnihilation);
  pManager->RemoveProcess(&thePositronIonisation);
  pManager->RemoveProcess(&thePositronMultipleScattering);
  if(synchOn) pManager->RemoveProcess(&thePositronSynch);
  }
  }
}

void G4EMBuilder::Build()
{
  G4ProcessManager * pManager = 0;
  wasActivated = true;
  
  pManager = G4Gamma::Gamma()->GetProcessManager();
  pManager->AddDiscreteProcess(&thePhotoEffect);
  pManager->AddDiscreteProcess(&theComptonEffect);
  pManager->AddDiscreteProcess(&thePairProduction);

  pManager = G4Electron::Electron()->GetProcessManager();
  pManager->AddDiscreteProcess(&theElectronBremsStrahlung);  
  pManager->AddProcess(&theElectronIonisation, ordInActive,2, 2);
  pManager->AddProcess(&theElectronMultipleScattering);
  if(synchOn) pManager->AddDiscreteProcess(&theElectronSynch);
  pManager->SetProcessOrdering(&theElectronMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&theElectronMultipleScattering, idxPostStep,  1);

  pManager = G4Positron::Positron()->GetProcessManager();
  pManager->AddDiscreteProcess(&thePositronBremsStrahlung);
  pManager->AddDiscreteProcess(&theAnnihilation);
  pManager->AddRestProcess(&theAnnihilation);
  pManager->AddProcess(&thePositronIonisation, ordInActive,2, 2);
  pManager->AddProcess(&thePositronMultipleScattering);
  if(synchOn) pManager->AddDiscreteProcess(&thePositronSynch);
  pManager->SetProcessOrdering(&thePositronMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&thePositronMultipleScattering, idxPostStep,  1);
  if (gammNucOn) theGNPhysics.Build();

}
// 2002 by J.P. Wellisch
