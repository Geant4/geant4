#ifndef G4NeutronHPMessenger_h
#define G4NeutronHPMessenger_h

#include "globals.hh"
#include "G4UImessenger.hh"

class G4NeutronHPManager;
class G4UIdirectory;
class G4UIcmdWithAString;

class G4NeutronHPMessenger: public G4UImessenger
{
   public:
      G4NeutronHPMessenger( G4NeutronHPManager* );
     ~G4NeutronHPMessenger();

      void SetNewValue(G4UIcommand*, G4String);

   private:
      G4NeutronHPManager* manager;

      G4UIdirectory* NeutronHPDir;
      G4UIcmdWithAString* PhotoEvaCmd;
      G4UIcmdWithAString* SkipMissingCmd;
      G4UIcmdWithAString* NeglectDopplerCmd;
      G4UIcmdWithAString* DoNotAdjustFSCmd;
      G4UIcmdWithAString* ProduceFissionFragementCmd;
      //G4UIcmdWithAString* AllowHeavyElementCmd;
/*
 * #setenv G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION 1
 * #setenv G4NEUTRONHP_SKIP_MISSING_ISOTOPES 1 
 * #setenv G4NEUTRONHP_NEGLECT_DOPPLER 1
 * #setenv G4NEUTRONHP_DO_NOT_ADJUST_FINAL_STATE 1
 * #setenv G4NEUTRONHP_PRODUCE_FISSION_FRAGMENTS 1
 *
*/
    
};

#endif
