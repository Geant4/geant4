// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserPhysicsListMessenger.hh,v 1.2 1999-04-14 10:45:55 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//---------------------------------------------------------------
//
//  G4UserPhysicsListMessenger.hh
//
//  Description:
//    This is a messenger class to interface to exchange information
//    between ParticleUserList and UI.
//
// ------------------------------------------------------
//  the List of Directory and Commands
// ------------------------------------------------------
//  /run/particle/   Paricle control commands.
//   Commands : 
//    SetCuts *  Set default cut value
//    DumpList * Dump List of particles in G4VUserPhysicsList.
//    dumpCutValue * Dump cut value information
//    Verbose * Set the Verbose level of G4VUserPhysicsList.
//    addProcessManager * add process manager
//    buildPhysicsTable * build physics table
// ------------------------------------------------------------
//	History
//        first version                   09 Jan. 1998 by H.Kurashige 
//        second version                  24 Jan. 1998 by H.Kurashige 
//        add buildPhysicsTable command   13 Apr. 1999 by H.Kurashige
// ------------------------------------------------------------

#ifndef G4UserPhysicsListMessenger_h
#define G4UserPhysicsListMessenger_h 1

class G4VUserPhysicsList;

class G4VUserPhysicsList;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString; 

#include "G4UImessenger.hh"
#include "globals.hh"

class G4UserPhysicsListMessenger: public G4UImessenger
{
  public:
    G4UserPhysicsListMessenger(G4VUserPhysicsList* pParticleList);
    virtual ~G4UserPhysicsListMessenger();
    
  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  protected:
    G4VUserPhysicsList* thePhysicsList;
    
  private: //commands
    G4UIdirectory *             theDirectory;
    G4UIcmdWithADoubleAndUnit * setCutCmd; 
    G4UIcmdWithAnInteger *      verboseCmd;
    G4UIcmdWithoutParameter *   dumpListCmd;
    G4UIcmdWithAString *        dumpCutValuesCmd;
    G4UIcmdWithAString *        addProcManCmd;
    G4UIcmdWithAString *        buildPTCmd;
};

#endif


