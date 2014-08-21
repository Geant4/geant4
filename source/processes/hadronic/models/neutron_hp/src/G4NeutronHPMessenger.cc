#include "G4NeutronHPMessenger.hh"
#include "G4NeutronHPManager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

G4NeutronHPMessenger::G4NeutronHPMessenger( G4NeutronHPManager* man )
:manager(man)
{
   NeutronHPDir = new G4UIdirectory( "/process/had/neutron_hp/" );
   NeutronHPDir->SetGuidance( "UI commands of NeutronHP" );

   PhotoEvaCmd = new G4UIcmdWithAString("/process/had/neutron_hp/use_photo_evaporation",this);
   PhotoEvaCmd->SetGuidance(" Force the use of the Photon Evaporation model, instead of the neutron capture final state data.");
   PhotoEvaCmd->SetParameterName("choice",false);
   PhotoEvaCmd->SetCandidates("true false");
   PhotoEvaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

   SkipMissingCmd = new G4UIcmdWithAString("/process/had/neutron_hp/skip_missing_isotopes",this);
   SkipMissingCmd->SetGuidance("Use only exact isotope data files, instead of allowing nearby isotope files to be used.");
   SkipMissingCmd->SetGuidance("In this case if the exact file is not available, the cross section will be set to zero.");
   SkipMissingCmd->SetParameterName("choice",false);
   SkipMissingCmd->SetCandidates("true false");
   SkipMissingCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

   NeglectDopplerCmd = new G4UIcmdWithAString("/process/had/neutron_hp/neglect_Doppler_broadening",this);
   NeglectDopplerCmd->SetGuidance("Switch off the Doppler broadening due to the thermal motion of the target nucleus.");
   NeglectDopplerCmd->SetGuidance("This option provides a significant CPU performance advantage.");
   NeglectDopplerCmd->SetParameterName("choice",false);
   NeglectDopplerCmd->SetCandidates("true false");
   NeglectDopplerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

   DoNotAdjustFSCmd = new G4UIcmdWithAString("/process/had/neutron_hp/do_not_adjust_final_state",this);
   DoNotAdjustFSCmd->SetGuidance("Disable to adjust final state for getting better conservation.");
   DoNotAdjustFSCmd->SetParameterName("choice",false);
   DoNotAdjustFSCmd->SetCandidates("true false");
   DoNotAdjustFSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

   ProduceFissionFragementCmd = new G4UIcmdWithAString("/process/had/neutron_hp/produce_fission_fragment",this);
   ProduceFissionFragementCmd->SetGuidance("Enable to generate fission fragments.");
   ProduceFissionFragementCmd->SetParameterName("choice",false);
   ProduceFissionFragementCmd->SetCandidates("true false");
   ProduceFissionFragementCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

G4NeutronHPMessenger::~G4NeutronHPMessenger()
{
   delete NeutronHPDir;
   delete PhotoEvaCmd;
   delete SkipMissingCmd;
   delete NeglectDopplerCmd;
   delete DoNotAdjustFSCmd;
   delete ProduceFissionFragementCmd;
}

void G4NeutronHPMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
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
}

