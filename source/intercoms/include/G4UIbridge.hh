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
// $Id: G4UIbridge.hh 66241 2012-12-13 18:34:42Z gunter $
//
// ====================================================================
//   G4UIbridge.hh
//
//  This is a concrete class of G4UIsession.
//
//  This class object is instantiated by G4UImanager at every time 
//  when "/control/execute macro_file" command is executed.
//  Also in the case of pure batch mode with a macro file, 
//  this class can be used as other ordinary G4UIsession 
//  concrete classes, i.e. SessionStart() is invoked in main().
// ====================================================================
#ifndef G4UIbridge_H
#define G4UIbridge_H 1

class G4UImanager;
#include "globals.hh"

// ====================================================================
//
// class definition
//
//  G4UIbridge:
//   To be used for MT mode. 
//   Register a particular thread-local G4UImanager with a UI command
//   directory name. When a UI command is issued in the master thread
//   that starts with this redistered directory name, it is immediately
//   forwarded to the registered G4UImanager. Such forwarded command 
//   is not processed in the master thread nor other worker thread.
//
// ====================================================================

class G4UIbridge 
{
 public:
  G4UIbridge(G4UImanager* localUI,G4String dir);
  ~G4UIbridge();

  G4int ApplyCommand(const G4String& aCmd);

 private:
  G4UImanager* localUImanager;
  G4String dirName;

 public:
  inline G4UImanager* LocalUI() const
  { return localUImanager; }
  inline G4String DirName() const
  { return dirName; }
  inline G4int DirLength() const
  { return dirName.length(); }
};

#endif
