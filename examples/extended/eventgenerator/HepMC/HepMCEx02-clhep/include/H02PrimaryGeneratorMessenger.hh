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
//   H02PrimaryGeneratorMessenger.hh
//   $Id: H02PrimaryGeneratorMessenger.hh,v 1.1 2002-11-19 10:35:49 murakami Exp $
//
// ====================================================================
#ifndef H02_PRIMARY_GENERATOR_MESSENGER_H
#define H02_PRIMARY_GENERATOR_MESSENGER_H

#include "globals.hh"
#include "G4UImessenger.hh"

class H02PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;

class H02PrimaryGeneratorMessenger : public G4UImessenger {
private:
  H02PrimaryGeneratorAction* primaryAction;

  G4UIdirectory* dir;
  G4UIcmdWithAString* select;

public:
  H02PrimaryGeneratorMessenger(H02PrimaryGeneratorAction* genaction);
  ~H02PrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand* command, G4String newValues);
  G4String GetCurrentValue(G4UIcommand* command);  
};


#endif
