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
// G4UIcontrolMessenger
//
// Class description:
//
// This class is a concrete class of G4UImessenger which defines
// commands affecting to the G4UImanager. Commands defined by
// this messenger are:
//   /control/
//   /control/macroPath
//   /control/execute
//   /control/loop
//   /control/foreach
//   /control/suppressAbortion
//   /control/verbose
//   /control/useDoublePrecision
//   /control/saveHistory
//   /control/stopSavingHistory
//   /control/alias
//   /control/unalias
//   /control/listAlias
//   /control/getEnv
//   /control/getVal
//   /control/echo
//   /control/shell
//   /control/manual
//   /control/createHTML
//   /control/maximumStoredHistory
//   /control/if
//   /control/doif
//   /control/add
//   /control/subtract
//   /control/multiply
//   /control/divide
//   /control/strif
//   /control/strdoif
//   /control/ifBatch
//   /control/ifInteractive
//   /control/doifBatch
//   /control/doifInteractive

// Author: Makoto Asai, SLAC - 2001
// --------------------------------------------------------------------
#ifndef G4UIcontrolMessenger_hh
#define G4UIcontrolMessenger_hh 1

#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcommand;

class G4UIcontrolMessenger : public G4UImessenger
{
  public:

    G4UIcontrolMessenger();
    ~G4UIcontrolMessenger() override;
    void SetNewValue(G4UIcommand* command, G4String newValue) override;
    G4String GetCurrentValue(G4UIcommand* command) override;

   private:

    G4UIdirectory* controlDirectory = nullptr;
    G4UIcmdWithAString* macroPathCommand = nullptr;
    G4UIcmdWithAString* ExecuteCommand = nullptr;
    G4UIcmdWithAnInteger* suppressAbortionCommand = nullptr;
    G4UIcmdWithAnInteger* verboseCommand = nullptr;
    G4UIcmdWithABool* doublePrecCommand = nullptr;
    G4UIcmdWithAString* historyCommand = nullptr;
    G4UIcmdWithoutParameter* stopStoreHistoryCommand = nullptr;
    G4UIcommand* aliasCommand = nullptr;
    G4UIcmdWithAString* unaliasCommand = nullptr;
    G4UIcmdWithoutParameter* listAliasCommand = nullptr;
    G4UIcmdWithAString* getEnvCmd = nullptr;
    G4UIcommand* getValCmd = nullptr;
    G4UIcmdWithAString* echoCmd = nullptr;
    G4UIcmdWithAString* shellCommand = nullptr;
    G4UIcommand* loopCommand = nullptr;
    G4UIcommand* foreachCommand = nullptr;
    G4UIcmdWithAString* ManualCommand = nullptr;
    G4UIcmdWithAString* HTMLCommand = nullptr;
    G4UIcmdWithAnInteger* maxStoredHistCommand = nullptr;
    G4UIcommand* ifCommand = nullptr;
    G4UIcommand* doifCommand = nullptr;
    G4UIcommand* addCommand = nullptr;
    G4UIcommand* subtractCommand = nullptr;
    G4UIcommand* multiplyCommand = nullptr;
    G4UIcommand* divideCommand = nullptr;
    G4UIcommand* remainderCommand = nullptr;
    G4UIcommand* strifCommand = nullptr;
    G4UIcommand* strdoifCommand = nullptr;
    G4UIcmdWithAString* ifBatchCommand = nullptr;
    G4UIcmdWithAString* ifInteractiveCommand = nullptr;
    G4UIcmdWithAString* doifBatchCommand = nullptr;
    G4UIcmdWithAString* doifInteractiveCommand = nullptr;
};

#endif
