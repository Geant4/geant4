// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProcessTableMessenger.hh,v 1.1 1999-01-07 16:13:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
//  G4ProcessTableMessenger.hh
//
//  Description:
//    This is a messenger class to interface to exchange information
//    between ProcessTable and UI.
//
//  /process/   Process Table control commands.
//   Commands : 
//     list *       :Dump process name registered.
//     verbose *    :Set Verbose Level for Process Table
//     activate *   :Activate process  
//     inactivate * :Inctivate process  
//     dump *       :Dump process information
//
//  History:
//    15 Aug. 1998, H. Kurashige  
//
//---------------------------------------------------------------

#ifndef G4ProcessTableMessenger_h
#define G4ProcessTableMessenger_h 1

class G4ProcessTable;

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger; 
class G4UIcmdWithAString;

#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4ProcessType.hh"

class G4ProcessTableMessenger: public G4UImessenger
{
  public:
    G4ProcessTableMessenger(G4ProcessTable* pTable);
    ~G4ProcessTableMessenger();

    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    G4String GetProcessTypeName(G4ProcessType aType) const;  
    G4int    GetProcessType(const G4String& aTypeName) const;
    static   G4int NumberOfProcessType;
    void     SetNumberOfProcessType();

  private:
    G4ProcessTableMessenger(const G4ProcessTableMessenger&){};
    G4ProcessTableMessenger(){};

  private:
    G4ProcessTable*  theProcessTable;

    G4UIdirectory *             thisDirectory;
    G4UIcmdWithAnInteger *      verboseCmd;
    G4UIcmdWithAString *        listCmd;
    G4UIcommand *               dumpCmd;
    G4UIcommand *               activateCmd;
    G4UIcommand *               inactivateCmd;

    G4String                    currentProcessTypeName;
    G4String                    currentProcessName;
    G4String                    currentParticleName; 
};

#endif




