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
