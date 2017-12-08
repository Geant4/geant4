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
#include "G4ParticleHPMessenger.hh"
#include "G4ParticleHPManager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

G4ParticleHPMessenger::G4ParticleHPMessenger( G4ParticleHPManager* man )
:manager(man)
{
   ParticleHPDir = new G4UIdirectory( "/process/had/particle_hp/" );
   ParticleHPDir->SetGuidance( "UI commands of ParticleHP" );

   PhotoEvaCmd = new G4UIcmdWithAString("/process/had/particle_hp/use_photo_evaporation",this);
   PhotoEvaCmd->SetGuidance(" Force the use of the Photon Evaporation model, instead of the neutron capture final state data.");
   PhotoEvaCmd->SetParameterName("choice",false);
   PhotoEvaCmd->SetCandidates("true false");
   PhotoEvaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

   SkipMissingCmd = new G4UIcmdWithAString("/process/had/particle_hp/skip_missing_isotopes",this);
   SkipMissingCmd->SetGuidance("Use only exact isotope data files, instead of allowing nearby isotope files to be used.");
   SkipMissingCmd->SetGuidance("In this case if the exact file is not available, the cross section will be set to zero.");
   SkipMissingCmd->SetParameterName("choice",false);
   SkipMissingCmd->SetCandidates("true false");
   SkipMissingCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

   NeglectDopplerCmd = new G4UIcmdWithAString("/process/had/particle_hp/neglect_Doppler_broadening",this);
   NeglectDopplerCmd->SetGuidance("Switch off the Doppler broadening due to the thermal motion of the target nucleus.");
   NeglectDopplerCmd->SetGuidance("This option provides a significant CPU performance advantage.");
   NeglectDopplerCmd->SetParameterName("choice",false);
   NeglectDopplerCmd->SetCandidates("true false");
   NeglectDopplerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

   DoNotAdjustFSCmd = new G4UIcmdWithAString("/process/had/particle_hp/do_not_adjust_final_state",this);
   DoNotAdjustFSCmd->SetGuidance("Disable to adjust final state for getting better conservation.");
   DoNotAdjustFSCmd->SetParameterName("choice",false);
   DoNotAdjustFSCmd->SetCandidates("true false");
   DoNotAdjustFSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

   ProduceFissionFragementCmd = new G4UIcmdWithAString("/process/had/particle_hp/produce_fission_fragment",this);
   ProduceFissionFragementCmd->SetGuidance("Enable to generate fission fragments.");
   ProduceFissionFragementCmd->SetParameterName("choice",false);
   ProduceFissionFragementCmd->SetCandidates("true false");
   ProduceFissionFragementCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

   NRESP71Cmd = new G4UIcmdWithAString("/process/had/particle_hp/use_NRESP71_model",this);
   NRESP71Cmd->SetGuidance("Enable to use NRESP71 model for n on C reaction");
   NRESP71Cmd->SetParameterName("choice",false);
   NRESP71Cmd->SetCandidates("true false");
   NRESP71Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

   VerboseCmd = new G4UIcmdWithAnInteger("/process/had/particle_hp/verbose",this);
   VerboseCmd->SetGuidance("Set Verbose level of ParticleHP package");
   VerboseCmd->SetParameterName("verbose_level",true);
   VerboseCmd->SetDefaultValue(1);
   VerboseCmd->SetRange("verbose_level >=0");
   VerboseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

G4ParticleHPMessenger::~G4ParticleHPMessenger()
{
   delete ParticleHPDir;
   delete PhotoEvaCmd;
   delete SkipMissingCmd;
   delete NeglectDopplerCmd;
   delete DoNotAdjustFSCmd;
   delete ProduceFissionFragementCmd;
   delete VerboseCmd;
}

void G4ParticleHPMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
   G4bool bValue=false;
   if ( newValue == "true" ) bValue=true;

   if ( command == PhotoEvaCmd ) { 
      manager->SetUseOnlyPhotoEvaporation( bValue ); 
   }
   if ( command == SkipMissingCmd) { 
      manager->SetSkipMissingIsotopes( bValue ); 
   }
   if ( command == NeglectDopplerCmd ) { 
      manager->SetNeglectDoppler( bValue ); 
   }
   if ( command == DoNotAdjustFSCmd ) { 
      manager->SetDoNotAdjustFinalState( bValue ); 
   }
   if ( command == ProduceFissionFragementCmd ) { 
      manager->SetProduceFissionFragments( bValue ); 
   }
   if ( command == NRESP71Cmd ) { 
      manager->SetUseNRESP71Model( bValue ); 
   }
   if ( command == VerboseCmd ) {
      manager->SetVerboseLevel( VerboseCmd->ConvertToInt( newValue ) ); 
   }
}

