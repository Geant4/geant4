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
//   ExN04PrimaryGeneratorMessenger.hh
//   $Id: ExN04PrimaryGeneratorMessenger.hh,v 1.1 2002-04-29 20:43:46 asaim Exp $
//
// ====================================================================
#ifndef EXN04_PRIMARY_GENERATOR_MESSENGER_H
#define EXN04_PRIMARY_GENERATOR_MESSENGER_H

#include "globals.hh"
#include "G4UImessenger.hh"

class ExN04PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;

class ExN04PrimaryGeneratorMessenger : public G4UImessenger {
private:
  ExN04PrimaryGeneratorAction* primaryAction;

  G4UIdirectory* dir;
  G4UIdirectory* mydetdir;
  G4UIcmdWithAString* select;

public:
  ExN04PrimaryGeneratorMessenger(ExN04PrimaryGeneratorAction* genaction);
  ~ExN04PrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand* command, G4String newValues);
  G4String GetCurrentValue(G4UIcommand* command);  
};


#endif
