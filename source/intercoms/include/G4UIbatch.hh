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
// $Id: G4UIbatch.hh 67965 2013-03-13 09:35:29Z gcosmo $
//
// ====================================================================
//   G4UIbatch.hh
//
//  This is a concrete class of G4UIsession.
//
//  This class object is instantiated by G4UImanager at every time 
//  when "/control/execute macro_file" command is executed.
//  Also in the case of pure batch mode with a macro file, 
//  this class can be used as other ordinary G4UIsession 
//  concrete classes, i.e. SessionStart() is invoked in main().
// ====================================================================
#ifndef G4UI_BATCH_H
#define G4UI_BATCH_H 1

#include "G4UIsession.hh"
#include <fstream>

// ====================================================================
//
// class definition
//
// ====================================================================

class G4UIbatch : public G4UIsession {
private:
  G4UIsession* previousSession;

  std::ifstream macroStream;
  G4bool isOpened;

  //static G4bool commandFailed;

  // get command from a batch script file
  G4String ReadCommand();
  G4int ExecCommand(const G4String& command);

public:
  G4UIbatch(const char* fileName, G4UIsession* prevSession=0);
  //  "prevSession" must be 0 if this class is constructed
  //  from main().

  ~G4UIbatch();
  
  G4UIsession* GetPreviousSession() const;

  virtual G4UIsession* SessionStart();
  virtual void PauseSessionStart(const G4String& Prompt);

};

// ============================================================================
// inlines

inline G4UIsession* G4UIbatch::GetPreviousSession() const
{
  return previousSession;
}

#endif
