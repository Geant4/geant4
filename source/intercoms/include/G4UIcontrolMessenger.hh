// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcontrolMessenger.hh,v 1.3 1999-12-15 14:50:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4UIcontrolMessenger_h
#define G4UIcontrolMessenger_h 1

#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;

// class description:
//  This class is a concrete class of G4UImessenger which defines
// commands affecting to the G4UImanager. Commands defined by
// this messenger are
//   /control/ 
//   /control/execute
//   /control/verbose
//   /control/saveHistory
//   /control/stopSavingHistory
//   /control/manual

class G4UIcontrolMessenger : public G4UImessenger 
{
  public:
      G4UIcontrolMessenger();
      ~G4UIcontrolMessenger();
      void SetNewValue(G4UIcommand * command,G4String newValue);
      G4String GetCurrentValue(G4UIcommand * command);

  private:
      G4UIdirectory * controlDirectory;
      G4UIcmdWithAString * ExecuteCommand;
      G4UIcmdWithAnInteger * verboseCommand;
      G4UIcmdWithAString * historyCommand;
      G4UIcmdWithoutParameter * stopStoreHistoryCommand;
      G4UIcmdWithAString * ManualCommand;
};

#endif

