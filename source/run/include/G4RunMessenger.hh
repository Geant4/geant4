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
// G4RunMessenger
//
// Class description:
//
// This is a messenger class for G4RunManager.
// Implemented commands are following;
//
//   beamOn           *  Start a Run.
//   verbose          *  Set the Verbose level of G4RunManager.
//   printProgress    *  Set the frequency of printing out progress of a run.
//   dumpRegion       *  Dump information of a region.
//   dumpCouples      *  Dump information of material-cuts-couples.
//   optimizeGeometry *  Set the optimization flag of closing geometry.
//   breakAtBeginOfEvent * Set a break point at the beginning of every event.
//   breakAtEndOfEvent   *   Set a break point at the end of every event.
//   abort            *  Abort current run processing.
//   Initialize       *  Initialise G4 kernel.
//   geometryModified *  Force geometry to be closed again.
//   physicsModified  *  Force cross-section tables to be calculated again
//                       (and rebuilding physics table will be invoked).
//   constructScoringWorlds * Construct scoring world(s) if defined.

// Original author: M.Asai, 1997
// --------------------------------------------------------------------
#ifndef G4RunMessenger_hh
#define G4RunMessenger_hh 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4RunManager;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
class G4UIcommand;
class G4MaterialScanner;

class G4RunMessenger : public G4UImessenger
{
  public:

    G4RunMessenger(G4RunManager* runMgr);
   ~G4RunMessenger();

    void SetNewValue(G4UIcommand* command, G4String newValues);
    G4String GetCurrentValue(G4UIcommand* command);

  private:

    G4RunManager* runManager = nullptr;
    G4String macroFileName = "***NULL***";  // internal use only!!!

    G4UIdirectory* runDirectory = nullptr;
    G4UIcommand* beamOnCmd = nullptr;
    G4UIcmdWithAnInteger* verboseCmd = nullptr;
    G4UIcmdWithAnInteger* printProgCmd = nullptr;
    G4UIcmdWithAnInteger* nThreadsCmd = nullptr;
    G4UIcmdWithoutParameter* maxThreadsCmd = nullptr;
    G4UIcmdWithAnInteger* pinAffinityCmd = nullptr;
    G4UIcommand* evModCmd = nullptr;
    G4UIcmdWithAString* dumpRegCmd = nullptr;
    G4UIcmdWithoutParameter* dumpCoupleCmd = nullptr;
    G4UIcmdWithABool* optCmd = nullptr;
    G4UIcmdWithABool* brkBoECmd = nullptr;
    G4UIcmdWithABool* brkEoECmd = nullptr;
    G4UIcmdWithABool* abortCmd = nullptr;
    G4UIcmdWithoutParameter* abortEventCmd = nullptr;
    G4UIcmdWithoutParameter* initCmd = nullptr;
    G4UIcmdWithoutParameter* geomCmd = nullptr;
    G4UIcmdWithABool* geomRebCmd = nullptr;
    G4UIcmdWithoutParameter* physCmd = nullptr;
    G4UIcmdWithAnInteger* randEvtCmd = nullptr;
    G4UIcommand* procUICmds = nullptr;

    G4UIdirectory* randomDirectory = nullptr;
    G4UIcmdWithAString* seedCmd = nullptr;
    G4UIcmdWithAString* randDirCmd = nullptr;
    G4UIcmdWithABool* savingFlagCmd = nullptr;
    G4UIcmdWithoutParameter* saveThisRunCmd = nullptr;
    G4UIcmdWithoutParameter* saveThisEventCmd = nullptr;
    G4UIcmdWithAString* restoreRandCmd = nullptr;
    G4UIcmdWithABool* saveEachEventCmd = nullptr;
    G4UIcmdWithABool* restoreRandCmdMT = nullptr;

    G4UIcmdWithoutParameter* constScoreCmd = nullptr;

    G4MaterialScanner* materialScanner = nullptr;
};

#endif
