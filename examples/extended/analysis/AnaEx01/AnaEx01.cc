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
// $Id: AnaEx01.cc,v 1.8 2001-11-16 14:35:10 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - AnaEx01
//
// --------------------------------------------------------------
// Comments
//   Example of histogram and tuple manipulations using an AIDA compliant 
//  system. All analysis manipulations (hooking an AIDA implementation,
//  histo booking, filling, etc...) are concentrated in one 
//  class : AnaEx01AnalysisManager.
//   See the README file within the same directory to have more infos.
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "Randomize.hh"

#include "AnaEx01AnalysisManager.hh"

#include "AnaEx01DetectorConstruction.hh"
#include "AnaEx01PhysicsList.hh"
#include "AnaEx01PrimaryGeneratorAction.hh"
#include "AnaEx01RunAction.hh"
#include "AnaEx01EventAction.hh"
#include "AnaEx01SteppingAction.hh"
#include "AnaEx01SteppingVerbose.hh"

int main(int argc,char** argv) {

  // choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new AnaEx01SteppingVerbose);
     
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  AnaEx01DetectorConstruction* detector = new AnaEx01DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new AnaEx01PhysicsList);

  runManager->SetUserAction(new AnaEx01PrimaryGeneratorAction(detector));

#ifdef G4ANALYSIS_USE
  AnaEx01AnalysisManager* analysisManager = new AnaEx01AnalysisManager();
  runManager->SetUserAction(new AnaEx01RunAction(analysisManager));
  runManager->SetUserAction(new AnaEx01EventAction(analysisManager));
  runManager->SetUserAction(new AnaEx01SteppingAction(analysisManager));
#else
  runManager->SetUserAction(new AnaEx01RunAction());
  runManager->SetUserAction(new AnaEx01EventAction());
  runManager->SetUserAction(new AnaEx01SteppingAction());
#endif

  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  // Batch mode
  UI->ApplyCommand("/control/execute run.mac");

  // job termination
#ifdef G4ANALYSIS_USE
  delete analysisManager;
#endif
  delete runManager;

  return 0;
}

