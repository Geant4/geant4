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
//
//---------------------------------------------------------------------------
//
// ClassName:   G4BertiniElectroNuclearBuilder
//
// Author: 2002 H.P. Wellisch
//
// Modified: 
// 25.04.2006 V.Ivanchenko fix destructor 
//
//----------------------------------------------------------------------------
//

#include "G4BertiniElectroNuclearBuilder.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"
#include "G4GammaGeneralProcess.hh"
#include "G4LossTableManager.hh"

#include "G4HadronicParameters.hh"

G4BertiniElectroNuclearBuilder::G4BertiniElectroNuclearBuilder(G4bool eNucl) : 
  thePhotoNuclearProcess(nullptr), theElectronNuclearProcess(nullptr), 
  thePositronNuclearProcess(nullptr), theElectroReaction(nullptr), 
  theGammaReaction(nullptr), theModel(nullptr), theCascade(nullptr), 
  theStringModel(nullptr), theFragmentation(nullptr), theStringDecay(nullptr), 
  wasActivated(false), eActivated(eNucl)
{}

G4BertiniElectroNuclearBuilder::~G4BertiniElectroNuclearBuilder() 
{
  if(wasActivated) {
    delete theFragmentation;
    delete theStringDecay;
  }
}

void G4BertiniElectroNuclearBuilder::Build()
{
  if(wasActivated) return;
  wasActivated=true;
  
  thePhotoNuclearProcess = new G4PhotoNuclearProcess;
  if(eActivated) {
    theElectronNuclearProcess = new G4ElectronNuclearProcess;
    thePositronNuclearProcess = new G4PositronNuclearProcess;
    theElectroReaction = new G4ElectroVDNuclearModel;
  }
  theGammaReaction = new G4CascadeInterface;

  theModel = new G4TheoFSGenerator;

  theStringModel = new G4QGSModel< G4GammaParticipants >;
  theStringDecay = new G4ExcitedStringDecay(theFragmentation=new G4QGSMFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);

  theCascade = new G4GeneratorPrecompoundInterface;

  theModel->SetTransport(theCascade);
  theModel->SetHighEnergyGenerator(theStringModel);

  G4ProcessManager * aProcMan = nullptr;

  theGammaReaction->SetMaxEnergy(3.5*GeV);
  thePhotoNuclearProcess->RegisterMe(theGammaReaction);
  theModel->SetMinEnergy(3.*GeV);
  theModel->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  thePhotoNuclearProcess->RegisterMe(theModel);
  
  G4GammaGeneralProcess* sp = 
    (G4GammaGeneralProcess*)G4LossTableManager::Instance()->GetGammaGeneralProcess();
  if(sp) {
    sp->AddHadProcess(thePhotoNuclearProcess);
  } else {
    aProcMan = G4Gamma::Gamma()->GetProcessManager();
    aProcMan->AddDiscreteProcess(thePhotoNuclearProcess);
  }

  if(eActivated) {  
    aProcMan = G4Electron::Electron()->GetProcessManager();
    theElectronNuclearProcess->RegisterMe(theElectroReaction);
    aProcMan->AddDiscreteProcess(theElectronNuclearProcess);
  
    aProcMan = G4Positron::Positron()->GetProcessManager();
    thePositronNuclearProcess->RegisterMe(theElectroReaction);
    aProcMan->AddDiscreteProcess(thePositronNuclearProcess);
  }
}

