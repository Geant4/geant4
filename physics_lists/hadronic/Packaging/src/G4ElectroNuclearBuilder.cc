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
#include "G4ElectroNuclearBuilder.hh"

#include "globals.hh"
#include "G4ios.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"

G4ElectroNuclearBuilder::G4ElectroNuclearBuilder() : wasActivated(false)
{
}

G4ElectroNuclearBuilder::~G4ElectroNuclearBuilder() 
{
  delete theFragmentation;
  delete theStringDecay;
  delete theStringModel;
  delete thePhotoNuclearProcess; 
  delete theElectronNuclearProcess;
  delete thePositronNuclearProcess;
  delete theElectroReaction;
  delete theGammaReaction;
  delete theModel;
  delete theCascade;
  if(wasActivated)
  {
    G4ProcessManager * pManager = 0;
    pManager = G4Gamma::Gamma()->GetProcessManager();
    if(pManager) pManager->RemoveProcess(thePhotoNuclearProcess);
    pManager = G4Electron::Electron()->GetProcessManager();
    if(pManager) pManager->RemoveProcess(theElectronNuclearProcess);
    pManager = G4Positron::Positron()->GetProcessManager();
    if(pManager) pManager->RemoveProcess(thePositronNuclearProcess);
  }
}

void G4ElectroNuclearBuilder::Build()
{
  wasActivated=true;
  
  thePhotoNuclearProcess = new G4PhotoNuclearProcess;
  theElectronNuclearProcess = new G4ElectronNuclearProcess;
  thePositronNuclearProcess = new G4PositronNuclearProcess;
  theElectroReaction = new G4ElectroNuclearReaction;
  theGammaReaction = new G4GammaNuclearReaction;

  theModel = new G4TheoFSGenerator;

  theStringModel = new G4QGSModel< G4GammaParticipants >;
  theStringDecay = new G4ExcitedStringDecay(theFragmentation=new G4QGSMFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);

  theCascade = new G4GeneratorPrecompoundInterface;

  theModel->SetTransport(theCascade);
  theModel->SetHighEnergyGenerator(theStringModel);

  G4ProcessManager * aProcMan = 0;
  
  aProcMan = G4Gamma::Gamma()->GetProcessManager();
  theGammaReaction->SetMaxEnergy(3.5*GeV);
  thePhotoNuclearProcess->RegisterMe(theGammaReaction);
  theModel->SetMinEnergy(3.*GeV);
  theModel->SetMaxEnergy(100*TeV);
  thePhotoNuclearProcess->RegisterMe(theModel);
  aProcMan->AddDiscreteProcess(thePhotoNuclearProcess);
  
  aProcMan = G4Electron::Electron()->GetProcessManager();
  theElectronNuclearProcess->RegisterMe(theElectroReaction);
  aProcMan->AddDiscreteProcess(theElectronNuclearProcess);
  
  aProcMan = G4Positron::Positron()->GetProcessManager();
  thePositronNuclearProcess->RegisterMe(theElectroReaction);
  aProcMan->AddDiscreteProcess(thePositronNuclearProcess);
}



// 2002 by J.P. Wellisch
