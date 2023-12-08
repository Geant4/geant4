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

#include "G4ElectroVDNuclearModel.hh"
#include "G4ElectronNuclearProcess.hh"
#include "G4PositronNuclearProcess.hh"
#include "G4PhotoNuclearCrossSection.hh"

#include "G4ios.hh"

#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSModel.hh"
#include "G4GammaParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"


#include "G4HadronicParameters.hh"

G4BertiniElectroNuclearBuilder::G4BertiniElectroNuclearBuilder(G4bool eNucl)
  : eActivated(eNucl)
{}

void G4BertiniElectroNuclearBuilder::Build()
{
  // gamma
  thePhotoNuclearProcess = new G4HadronInelasticProcess( "photonNuclear", G4Gamma::Definition() );
  thePhotoNuclearProcess->AddDataSet( new G4PhotoNuclearCrossSection() );
  theGammaReaction = new G4CascadeInterface();

  auto theModel = new G4TheoFSGenerator;

  auto theStringModel = new G4QGSModel< G4GammaParticipants >;
  auto theStringDecay = new G4ExcitedStringDecay( new G4QGSMFragmentation() );
  theStringModel->SetFragmentationModel(theStringDecay);

  auto theCascade = new G4GeneratorPrecompoundInterface();

  theModel->SetTransport(theCascade);
  theModel->SetHighEnergyGenerator(theStringModel);

  G4ProcessManager* aProcMan = nullptr;

  theGammaReaction->SetMaxEnergy(3.5*GeV);
  thePhotoNuclearProcess->RegisterMe(theGammaReaction);
  theModel->SetMinEnergy(3.*GeV);
  theModel->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  thePhotoNuclearProcess->RegisterMe(theModel);
  
  G4GammaGeneralProcess* sp = 
    dynamic_cast<G4GammaGeneralProcess*>(G4LossTableManager::Instance()->GetGammaGeneralProcess());
  if ( nullptr != sp ) {
    sp->AddHadProcess(thePhotoNuclearProcess);
  } else {
    aProcMan = G4Gamma::Gamma()->GetProcessManager();
    aProcMan->AddDiscreteProcess(thePhotoNuclearProcess);
  }

  // e+, e-
  if (eActivated) {  

    auto theElectronNuclearProcess = new G4ElectronNuclearProcess();
    auto thePositronNuclearProcess = new G4PositronNuclearProcess();
    auto theElectroReaction = new G4ElectroVDNuclearModel();

    aProcMan = G4Electron::Electron()->GetProcessManager();
    theElectronNuclearProcess->RegisterMe(theElectroReaction);
    aProcMan->AddDiscreteProcess(theElectronNuclearProcess);
  
    aProcMan = G4Positron::Positron()->GetProcessManager();
    thePositronNuclearProcess->RegisterMe(theElectroReaction);
    aProcMan->AddDiscreteProcess(thePositronNuclearProcess);
  }
}

