#include "G4ElectroNuclearBuilder.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/iomanip"   

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"

G4ElectroNuclearBuilder::G4ElectroNuclearBuilder() 
{
  theElectroReaction = new G4ElectroNuclearReaction;
  theGammaReaction = new G4GammaNuclearReaction;
  theModel = new G4TheoFSGenerator;
  theCascade = new G4StringChipsParticleLevelInterface;
  theModel->SetTransport(theCascade);
  theModel->SetHighEnergyGenerator(&theStringModel);
  theStringDecay = new G4ExcitedStringDecay(&theFragmentation);
  theStringModel.SetFragmentationModel(theStringDecay);
}

G4ElectroNuclearBuilder::~G4ElectroNuclearBuilder() {delete theStringDecay;}

void G4ElectroNuclearBuilder::Build()
{
  G4ProcessManager * aProcMan = 0;
  
  aProcMan = G4Gamma::Gamma()->GetProcessManager();
  theGammaReaction->SetMaxEnergy(3.5*GeV);
  thePhotoNuclearProcess.RegisterMe(theGammaReaction);
  theModel->SetMinEnergy(3.*GeV);
  theModel->SetMaxEnergy(100*TeV);
  thePhotoNuclearProcess.RegisterMe(theModel);
  aProcMan->AddDiscreteProcess(&thePhotoNuclearProcess);
  
  aProcMan = G4Electron::Electron()->GetProcessManager();
  theElectronNuclearProcess.RegisterMe(theElectroReaction);
  aProcMan->AddDiscreteProcess(&theElectronNuclearProcess);
  
  aProcMan = G4Positron::Positron()->GetProcessManager();
  thePositronNuclearProcess.RegisterMe(theElectroReaction);
  aProcMan->AddDiscreteProcess(&thePositronNuclearProcess);
}



// 2002 by J.P. Wellisch
