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




