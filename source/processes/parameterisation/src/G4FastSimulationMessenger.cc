// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FastSimulationMessenger.cc,v 1.2 1999-04-14 14:25:35 mora Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

  fListEnvelopesCmd = 
    new G4UIcmdWithAString("/param/listEnvelopes", this);
  fListEnvelopesCmd->SetParameterName("ParticleName",true);
  fListEnvelopesCmd->SetDefaultValue("all");
  fListEnvelopesCmd->SetGuidance("List all the envelope names for a given Particle");
  fListEnvelopesCmd->SetGuidance("(or for all particles if without parameters).");
  fListEnvelopesCmd->AvailableForStates(PreInit,Idle);
  
  fListModelsCmd = 
    new G4UIcmdWithAString("/param/listModels", this);
  fListModelsCmd->SetParameterName("EnvelopeName",true);
  fListModelsCmd->SetDefaultValue("all");
  fListModelsCmd->SetGuidance("List all the Model names for a given Envelope");
  fListModelsCmd->SetGuidance("(or for all envelopes if without parameters).");
  fListModelsCmd->AvailableForStates(PreInit,Idle);
  
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
  if( command == fListEnvelopesCmd)
    if(newValue == "all") 
      fGlobalFastSimulationManager->ListEnvelopes();
    else 
      fGlobalFastSimulationManager->
	ListEnvelopes(G4ParticleTable::GetParticleTable()->
		      FindParticle(newValue));
  if( command == fListModelsCmd)
    fGlobalFastSimulationManager->ListEnvelopes(newValue, MODELS); 
  if( command == fListIsApplicableCmd)
    fGlobalFastSimulationManager->ListEnvelopes(newValue, ISAPPLICABLE);
  if( command == fActivateModel)
    fGlobalFastSimulationManager->ActivateFastSimulationModel(newValue);
  if( command == fInActivateModel)
    fGlobalFastSimulationManager->InActivateFastSimulationModel(newValue);
}
