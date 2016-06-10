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
// $Id: G4FastSimulationMessenger.cc 68056 2013-03-13 14:44:48Z gcosmo $
//

#include "G4FastSimulationMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "G4ios.hh"

G4FastSimulationMessenger::
G4FastSimulationMessenger(G4GlobalFastSimulationManager* theGFSM)
  : fGlobalFastSimulationManager(theGFSM)
{
  fFSDirectory = new G4UIdirectory("/param/");
  fFSDirectory->SetGuidance("Fast Simulation print/control commands.");

  fShowSetupCmd =
    new G4UIcmdWithoutParameter("/param/showSetup", this);
  fShowSetupCmd->SetGuidance("Show fast simulation setup:");
  fShowSetupCmd->SetGuidance("    - for each world region:");
  fShowSetupCmd->SetGuidance("        1) fast simulation manager process attached;");
  fShowSetupCmd->SetGuidance("               - and to which particles the process is attached to;");
  fShowSetupCmd->SetGuidance("        2) region hierarchy;");
  fShowSetupCmd->SetGuidance("               - with for each the fast simulation models attached;");
  fShowSetupCmd->AvailableForStates(G4State_Idle, G4State_GeomClosed);

  fListEnvelopesCmd = 
    new G4UIcmdWithAString("/param/listEnvelopes", this);
  fListEnvelopesCmd->SetParameterName("ParticleName",true);
  fListEnvelopesCmd->SetDefaultValue("all");
  fListEnvelopesCmd->SetGuidance("List all the envelope names for a given Particle");
  fListEnvelopesCmd->SetGuidance("(or for all particles if without parameters).");
  fListEnvelopesCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fListModelsCmd = 
    new G4UIcmdWithAString("/param/listModels", this);
  fListModelsCmd->SetParameterName("EnvelopeName",true);
  fListModelsCmd->SetDefaultValue("all");
  fListModelsCmd->SetGuidance("List all the Model names for a given Envelope");
  fListModelsCmd->SetGuidance("(or for all envelopes if without parameters).");
  fListModelsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fListIsApplicableCmd =
    new G4UIcmdWithAString("/param/listIsApplicable", this);
  fListIsApplicableCmd->SetParameterName("ModelName",true);
  fListIsApplicableCmd->SetDefaultValue("all");
  fListIsApplicableCmd->SetGuidance("List all the Particle names a given Model is applicable");
  fListIsApplicableCmd->SetGuidance("(or for all Models if without parameters).");

  fActivateModel =
    new G4UIcmdWithAString("/param/ActivateModel", this);
  fActivateModel->SetParameterName("ModelName",false);
  fActivateModel->SetGuidance("Activate a given Model.");

  fInActivateModel =
    new G4UIcmdWithAString("/param/InActivateModel", this);
  fInActivateModel->SetParameterName("ModelName",false);
  fInActivateModel->SetGuidance("InActivate a given Model.");
}

G4FastSimulationMessenger::~G4FastSimulationMessenger()
{
  delete fShowSetupCmd;
  fShowSetupCmd = 0;
  delete fListIsApplicableCmd;
  fListIsApplicableCmd = 0;
  delete fActivateModel;
  fActivateModel = 0;
  delete fInActivateModel;
  fInActivateModel = 0;
  delete fListModelsCmd;
  fListModelsCmd = 0;
  delete fListEnvelopesCmd;
  fListEnvelopesCmd = 0;
  delete fFSDirectory;
  fFSDirectory = 0;
}

void G4FastSimulationMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if (command == fShowSetupCmd)
    fGlobalFastSimulationManager->ShowSetup();
  if( command == fListEnvelopesCmd)
  {
    if(newValue == "all") 
      fGlobalFastSimulationManager->ListEnvelopes();
    else 
      fGlobalFastSimulationManager->
	ListEnvelopes(G4ParticleTable::GetParticleTable()->
		      FindParticle(newValue));
  }
  if( command == fListModelsCmd)
    fGlobalFastSimulationManager->ListEnvelopes(newValue, MODELS); 
  if( command == fListIsApplicableCmd)
    fGlobalFastSimulationManager->ListEnvelopes(newValue, ISAPPLICABLE);
  if( command == fActivateModel)
    fGlobalFastSimulationManager->ActivateFastSimulationModel(newValue);
  if( command == fInActivateModel)
    fGlobalFastSimulationManager->InActivateFastSimulationModel(newValue);
}
