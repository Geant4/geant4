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
// $Id: G4RunMessenger.hh,v 1.18 2007-11-13 15:48:44 asaim Exp $
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
//    dumpRegion *          Dump information of a region.
//    dumpCouples *         Dump information of material-cuts-couples.
//    optimizeGeometry *    Set the optimization flag of closing geometry.
//    breakAtBeginOfEvent * Set a break point at the begining of every event.
//    breakAtEndOfEvent *   Set a break point at the end of every event.
//    abort *               Abort current run processing.
//    Initialize *          Initialize G4 kernel.
//    geometryModified *    Force geometry to be closed again.
//    physicsModified *     Force cross-section tables to be calculated again.
//           (and rebuilding physics table will be invoked)
//    constructScoringWorlds * Constrct scoring world(s) if defined
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
class G4MaterialScanner;

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
    G4UIcmdWithAString *        dumpRegCmd;
    G4UIcmdWithoutParameter *   dumpCoupleCmd;
    G4UIcmdWithABool *          optCmd;
    G4UIcmdWithABool *          brkBoECmd;
    G4UIcmdWithABool *          brkEoECmd;
    G4UIcmdWithABool *          abortCmd;
    G4UIcmdWithoutParameter *   abortEventCmd;
    G4UIcmdWithoutParameter *   initCmd;
    G4UIcmdWithoutParameter *   geomCmd;
    G4UIcmdWithoutParameter *   physCmd;
    G4UIcmdWithoutParameter *   cutCmd;
    G4UIcmdWithAnInteger *      randEvtCmd;

    G4UIdirectory *             randomDirectory;
    G4UIcmdWithAString *        seedCmd;
    G4UIcmdWithAString *        randDirCmd;
    G4UIcmdWithABool *          savingFlagCmd;
    G4UIcmdWithoutParameter *   saveThisRunCmd;
    G4UIcmdWithoutParameter *   saveThisEventCmd;
    G4UIcmdWithAString *        restoreRandCmd;
    
    G4UIcmdWithAString *        randDirOld;
    G4UIcmdWithAnInteger *      storeRandOld;
    G4UIcmdWithAString *        restoreRandOld;  

    G4UIcmdWithoutParameter *   constScoreCmd;

    G4MaterialScanner *         materialScanner;
};

#endif


