//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************

// ====================================================================
//
//   XXPrimaryGeneratorMessenger.hh
//   $Id: XXPrimaryGeneratorMessenger.hh,v 1.1 2002-04-29 20:44:45 asaim Exp $
//
// ====================================================================
#ifndef XX_PRIMARY_GENERATOR_MESSENGER_H
#define XX_PRIMARY_GENERATOR_MESSENGER_H

#include "globals.hh"
#include "G4UImessenger.hh"

class XXPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;

class XXPrimaryGeneratorMessenger : public G4UImessenger {
private:
  XXPrimaryGeneratorAction* primaryAction;

  G4UIdirectory* dir;
  G4UIcmdWithAString* select;

public:
  XXPrimaryGeneratorMessenger(XXPrimaryGeneratorAction* genaction);
  ~XXPrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand* command, G4String newValues);
  G4String GetCurrentValue(G4UIcommand* command);  
};


#endif
