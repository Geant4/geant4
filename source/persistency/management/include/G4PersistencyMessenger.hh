// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistencyMessenger.hh,v 1.4 1999/11/26 11:19:39 morita Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// 
//	GEANT 4 class header file 

// class description:
//
//      This is a messenger class for G4PersistencyManager.
//      Implemented commands are following;
//
//  Commands : 
//    /db/verbose *         Set the Verbose level of G4PersistencyManager.
//    /db/run *             Set the name of the Run Database
//    /db/event *           Set the name of the Event Database
//    /db/geometry *        Set the name of the Geometry Database
// 

#ifndef G4PersistencyMessenger_h
#define G4PersistencyMessenger_h 1

class G4PersistencyManager;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
class G4UIcommand;

#include "G4UImessenger.hh"
#include "globals.hh"

class G4PersistencyMessenger: public G4UImessenger
{
  public:
    G4PersistencyMessenger(G4PersistencyManager* persistencyMgr);
    ~G4PersistencyMessenger();

  public:
    void SetNewValue(G4UIcommand* command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand* command);

  private:
    G4PersistencyManager* persistencyManager;

  private: //commands
    G4UIdirectory*             persistencyDirectory;
    G4UIcmdWithAnInteger*      verboseCmd;
    G4UIcmdWithAString*        runDbCmd;   
    G4UIcmdWithAString*        eventDbCmd;   
    G4UIcmdWithAString*        geomDbCmd;   
};

#endif

