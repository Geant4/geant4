//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4RunMessenger.hh,v 1.9 2002-08-08 17:29:26 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//	GEANT 4 class header file 

// class description:
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
//    storeRandomNumberStatus *   Set the flag for storing random number status
//    restoreRandomNumberStatus * Restore the stored random number status
// 

#ifndef G4RunMessenger_h
#define G4RunMessenger_h 1

class G4RunManager;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
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
    G4UIcmdWithABool *          abortCmd;
    G4UIcmdWithoutParameter *   abortEventCmd;
    G4UIcmdWithoutParameter *   initCmd;
    G4UIcmdWithoutParameter *   geomCmd;
    G4UIcmdWithoutParameter *   cutCmd;
    
    G4UIdirectory *             randomDirectory;
    G4UIcmdWithAString *        randDirCmd;
    G4UIcmdWithABool *          savingFlagCmd;
    G4UIcmdWithoutParameter *   saveThisRunCmd;
    G4UIcmdWithoutParameter *   saveThisEventCmd;
    G4UIcmdWithAString *        restoreRandCmd;
    
    G4UIcmdWithAString *        randDirOld;
    G4UIcmdWithAnInteger *      storeRandOld;
    G4UIcmdWithAString *        restoreRandOld;  
};

#endif


