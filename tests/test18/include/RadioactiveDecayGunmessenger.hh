#ifndef RadioactiveDecayGunmessenger_h
#define RadioactiveDecayGunmessenger_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              RadioactiveDecayGunmessenger.hh
//
// Version:             0.b.3
// Date:                29/02/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "G4UImessenger.hh"

#include "UIcmdWithNucleusAndUnit.hh"
#include "globals.hh"

#include "RadioactiveDecayGun.hh"

class RadioactiveDecayGun;
////////////////////////////////////////////////////////////////////////////////
//
class RadioactiveDecayGunmessenger: public G4UImessenger
{
  // class description
  //   The RadioactiveDecayGunmessenger is instatiated by the
  //   RadioactiveDecayGun and introduces into the UI menu commands options
  //   to allow user to specify the isotope which will be treated by
  //   G4RadioactiveDecay
  // class description - end

public: // with description
  RadioactiveDecayGunmessenger (RadioactiveDecayGun *theRadioactiveDecayGun1);
  //    Constructor introduces commands into the UI menu to control
  //    RadioactiveDecayGun.  theG4RadioactiveDecayGun1 is used to identify
  //    to this class (when instatiated) the associated RadioactiveDecayGun,
  //    whose parameters are going to be changed as a result of the UI commands.
  //
  ~RadioactiveDecayGunmessenger ();
  //    Destructor deletes G4UIdirectory and G4UIcommand objects.
  //
  void SetNewValue (G4UIcommand *command, G4String newValues);
  //    Identifies the command which has been invoked by the user, extracts the
  //    parameters associated with that command (held in newValues, and uses
  //    these values with the appropriate member function of
  //    RadioactiveDecayGun.
  //
private:
  RadioactiveDecayGun          *theRadioactiveDecayGun;
  UIcmdWithNucleusAndUnit      *ionCmd;
};
////////////////////////////////////////////////////////////////////////////////
#endif




