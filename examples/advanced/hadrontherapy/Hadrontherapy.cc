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
// $Id: Hadrontherapy.cc Main of the Hadrontherapy example; 
// Last modified: G.A.P.Cirrone 
// 
// See more at: http://workgroup.lngs.infn.it/geant4lns/
//
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova, Genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "HadrontherapyEventAction.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyDetectorSD.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapyRunAction.hh"
#include "HadrontherapyMatrix.hh"
#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UImessenger.hh"
#include "globals.hh"
#include "HadrontherapySteppingAction.hh"

#include "HadrontherapyAnalysisManager.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

//////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc ,char ** argv)
{
  // Set the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());


  G4RunManager* runManager = new G4RunManager;

  // Initialize the geometry
  runManager -> SetUserInitialization(new HadrontherapyDetectorConstruction());

  // Initialize the physics 
  runManager -> SetUserInitialization(new HadrontherapyPhysicsList());

  // Initialize the primary particles
  HadrontherapyPrimaryGeneratorAction *pPrimaryGenerator = new HadrontherapyPrimaryGeneratorAction();
  runManager -> SetUserAction(pPrimaryGenerator);

  // Initialize matrix 
  HadrontherapyMatrix* matrix = HadrontherapyMatrix::getInstance();
  matrix -> Initialize();

  // Optional UserActions: run, event, stepping
  HadrontherapyRunAction* pRunAction = new HadrontherapyRunAction();
  runManager -> SetUserAction(pRunAction);

  HadrontherapyEventAction* pEventAction = new HadrontherapyEventAction(matrix);
  runManager -> SetUserAction(pEventAction);

  HadrontherapySteppingAction* steppingAction = new HadrontherapySteppingAction(pRunAction); 
  runManager -> SetUserAction(steppingAction);    

#ifdef ANALYSIS_USE
  HadrontherapyAnalysisManager* analysis = 
    HadrontherapyAnalysisManager::getInstance();
  analysis -> book();
#endif

#ifdef G4VIS_USE
  // Visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager -> Initialize();
#endif 

  // Define UI session for interactive mode.
  G4UIsession* session = 0;
  if (argc == 1)   
    {
      session = new G4UIterminal();
    }

  // Get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  
  if (session)   // Define UI session for interactive mode.
    { 
      G4cout<<" UI session starts ..."<< G4endl;
      UI -> ApplyCommand("/control/execute defaultMacro.mac");    
      session -> SessionStart();
      delete session;
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI -> ApplyCommand(command + fileName);
    }

  matrix -> TotalEnergyDeposit();

#ifdef ANALYSIS_USE
  analysis -> finish();
#endif

  // Job termination
#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;

  return 0;
}
