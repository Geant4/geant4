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
// -------------------------------------------------------------------
// -------------------------------------------------------------------

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;
class TrackingAction;

class RunAction : public G4UserRunAction
{
public:
  
  RunAction();
  ~RunAction();

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

private:

  /////////////////
  // Histogramming
  //
  void CreateHistogram();
  void WriteHistogram();

  /////////////////
  // Worker
  //
  void BeginMaster(const G4Run*);
  void EndMaster(const G4Run*);

  /////////////////
  // Worker
  //
  void InitializeWorker(const G4Run*);
  void BeginWorker(const G4Run*);
  void EndWorker(const G4Run*);

  /////////////////
  // Print Info
  //
  void PrintRunInfo(const G4Run* run);

  /////////////////
  // Attributes
  //
  TrackingAction* fpTrackingAction;
  G4String fFileName; 
  bool fInitialized;
  bool fDebug;
  
};
#endif
