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
//
// $Id: G4UIcontrolMessenger.hh 102561 2017-02-09 08:16:05Z gcosmo $
//

#ifndef G4UIcontrolMessenger_h
#define G4UIcontrolMessenger_h 1

#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcommand;

// class description:
//  This class is a concrete class of G4UImessenger which defines
// commands affecting to the G4UImanager. Commands defined by
// this messenger are
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

class G4UIcontrolMessenger : public G4UImessenger 
{
  public:
      G4UIcontrolMessenger();
      ~G4UIcontrolMessenger();
      void SetNewValue(G4UIcommand * command,G4String newValue);
      G4String GetCurrentValue(G4UIcommand * command);

  private:
      G4UIdirectory * controlDirectory;
      G4UIcmdWithAString * macroPathCommand;
      G4UIcmdWithAString * ExecuteCommand;
      G4UIcmdWithAnInteger * suppressAbortionCommand;
      G4UIcmdWithAnInteger * verboseCommand;
      G4UIcmdWithABool * doublePrecCommand;
      G4UIcmdWithAString * historyCommand;
      G4UIcmdWithoutParameter * stopStoreHistoryCommand;
      G4UIcommand * aliasCommand;
      G4UIcmdWithAString * unaliasCommand;
      G4UIcmdWithoutParameter * listAliasCommand;
      G4UIcmdWithAString * getEnvCmd;
      G4UIcommand * getValCmd;
      G4UIcmdWithAString * echoCmd;
      G4UIcmdWithAString * shellCommand;
      G4UIcommand * loopCommand;
      G4UIcommand * foreachCommand;
      G4UIcmdWithAString * ManualCommand;
      G4UIcmdWithAString * HTMLCommand;
      G4UIcmdWithAnInteger * maxStoredHistCommand;
      G4UIcommand * ifCommand;
      G4UIcommand * doifCommand;
      G4UIcommand * addCommand;
      G4UIcommand * subtractCommand;
      G4UIcommand * multiplyCommand;
      G4UIcommand * divideCommand;
      G4UIcommand * remainderCommand;
      G4UIcommand * strifCommand;
      G4UIcommand * strdoifCommand;
      G4UIcmdWithAString * ifBatchCommand;
      G4UIcmdWithAString * ifInteractiveCommand;
      G4UIcmdWithAString * doifBatchCommand;
      G4UIcmdWithAString * doifInteractiveCommand;
};

#endif

