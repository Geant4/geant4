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
//
//
//

#ifndef G4HumanPhantomPrimaryGeneratorMessenger_h
#define G4HumanPhantomPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4HumanPhantomPrimaryGeneratorAction;

class G4UIcmdWithAString;

class G4HumanPhantomPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
  G4HumanPhantomPrimaryGeneratorMessenger(G4HumanPhantomPrimaryGeneratorAction*);
  ~G4HumanPhantomPrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  G4HumanPhantomPrimaryGeneratorAction*    primary; 
  G4UIcmdWithAString*               beamCmd;
};

#endif

