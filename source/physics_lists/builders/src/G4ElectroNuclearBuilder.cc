//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4ElectroNuclearBuilder
//
// Author: 2002 H.P. Wellisch
//
// Modified: 
// 25.04.2006 V.Ivanchenko fix destructor 
//
//----------------------------------------------------------------------------
//

#include "G4ElectroNuclearBuilder.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"

G4ElectroNuclearBuilder::G4ElectroNuclearBuilder() : 
    thePhotoNuclearProcess(0), theElectronNuclearProcess(0), 
    thePositronNuclearProcess(0), theElectroReaction(0), 
    theGammaReaction(0), theModel(0), theCascade(0), 
    theStringModel(0), theFragmentation(0), theStringDecay(0), 
   wasActivated(false)
{
}

G4ElectroNuclearBuilder::~G4ElectroNuclearBuilder() 
{
  if(wasActivated) {
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
  }
}

void G4ElectroNuclearBuilder::Build()
{
  if(wasActivated) return;
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

