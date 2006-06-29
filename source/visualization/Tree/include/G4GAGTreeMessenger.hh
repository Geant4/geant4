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
