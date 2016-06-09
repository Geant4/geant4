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
// $Id: AnaEx01.cc,v 1.13 2005/12/06 11:07:17 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

// Geant4 :
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

// AIDA :
#ifdef G4ANALYSIS_USE
#include <AIDA/IAnalysisFactory.h>
#endif

// AnaEx01 :
#include "AnaEx01AnalysisManager.hh"
#include "AnaEx01DetectorConstruction.hh"
#include "AnaEx01PhysicsList.hh"
#include "AnaEx01PrimaryGeneratorAction.hh"
#include "AnaEx01RunAction.hh"
#include "AnaEx01EventAction.hh"
#include "AnaEx01SteppingAction.hh"
#include "AnaEx01SteppingVerbose.hh"

int main(int,char**) {

  // choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new AnaEx01SteppingVerbose);
     
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  AnaEx01DetectorConstruction* detector = new AnaEx01DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new AnaEx01PhysicsList);

  runManager->SetUserAction(new AnaEx01PrimaryGeneratorAction(detector));

  AnaEx01AnalysisManager* analysisManager = 0;
#ifdef G4ANALYSIS_USE
  AIDA::IAnalysisFactory* aida = AIDA_createAnalysisFactory();
  analysisManager = new AnaEx01AnalysisManager(aida);
#endif
  runManager->SetUserAction(new AnaEx01RunAction(analysisManager));
  runManager->SetUserAction(new AnaEx01EventAction(analysisManager));
  runManager->SetUserAction(new AnaEx01SteppingAction(analysisManager));

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

#ifdef G4ANALYSIS_USE
  delete aida;
#endif

  return 0;
}

