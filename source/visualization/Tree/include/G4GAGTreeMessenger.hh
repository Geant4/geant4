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
// Satoshi Tanaka 31th May 2001
//
// A messenger for G4GAGTree driver.


#ifndef G4GAGTREEMESSENGER_HH
#define G4GAGTREEMESSENGER_HH

#include "G4UImessenger.hh"

class G4UIcommand;
class G4UIcmdWithAnInteger;
class G4GAGTree;

class G4GAGTreeMessenger: public G4UImessenger {
public:
  G4GAGTreeMessenger(G4GAGTree*);
  virtual ~G4GAGTreeMessenger();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4GAGTree* fpGAGTree;
  G4UIcommand* fpDirectory;
  G4UIcmdWithAnInteger* fpCommandVerbose;
};

#endif
