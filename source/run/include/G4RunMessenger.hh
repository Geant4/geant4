// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RunMessenger.hh,v 1.1 1999-01-07 16:14:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//	GEANT 4 class header file 
//
//      This is a messenger class for G4RunManager.
//      Implemented commands are following;
//
//  Commands : 
//    BeamOn *              Start a Run.
//    Verbose *             Set the Verbose level of G4RunManager.
//    optimizeGeometry *    Set the optimization flag of closing geometry.
//    breakAtBeginOfEvent * Set a break point at the begining of every event.
//    breakAtEndOfEvent *   Set a break point at the end of every event.
//    abort *               Abort current run processing.
//    Initialize *          Initialize G4 kernel.
//    geometryModified *    Force geometry to be closed again.
//    cutoffModified *      Force closssection tables to be calculated again.
//           (and rebuilding physics table will be invoked)
// 
//	History
//        first version                              by Makoto Asai 
// ------------------------------------------------------------

#ifndef G4RunMessenger_h
#define G4RunMessenger_h 1

class G4RunManager;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
class G4UIcommand;

#include "G4UImessenger.hh"
#include "globals.hh"

class G4RunMessenger: public G4UImessenger
{
  public:
    G4RunMessenger(G4RunManager* runMgr);
    ~G4RunMessenger();
    
  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    G4RunManager * runManager;
    G4String macroFileName; // internal use only!!!
    
  private: //commands
    G4UIdirectory *             runDirectory;
    G4UIcommand *               beamOnCmd;
    G4UIcmdWithAnInteger *      verboseCmd;
    G4UIcmdWithABool *          optCmd;
    G4UIcmdWithABool *          brkBoECmd;
    G4UIcmdWithABool *          brkEoECmd;
    G4UIcmdWithoutParameter *   abortCmd;
    G4UIcmdWithoutParameter *   initCmd;
    G4UIcmdWithoutParameter *   geomCmd;
    G4UIcmdWithoutParameter *   cutCmd;
    
};

#endif


