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
#include "HadrontherapyInteractionParameters.hh"
#include "HadrontherapyAnalysisManager.hh"
#include "IAEADetectorConstruction.hh"
#include "HadrontherapyGeometryController.hh"
#include "HadrontherapyGeometryMessenger.hh"

#if defined(G4UI_USE_TCSH)
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#endif

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE_QT
#include "G4UIQt.hh"
#include "G4Qt.hh"
#endif

//////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc ,char ** argv)
{
  // Set the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());

  G4RunManager* runManager = new G4RunManager;
  //Initialize possible analysis needs, needs to come early in order to pick up metadata
#ifdef ANALYSIS_USE
  HadrontherapyAnalysisManager* analysis = HadrontherapyAnalysisManager::getInstance();
  analysis -> book();
#endif
  
  // Initialize the geometry user interface
  HadrontherapyGeometryController *geometryController = new HadrontherapyGeometryController();
  HadrontherapyGeometryMessenger *geometryMessenger = new HadrontherapyGeometryMessenger(geometryController);
  
  // Initialize the geometry
  HadrontherapyDetectorConstruction* pDetect = new HadrontherapyDetectorConstruction();
  runManager -> SetUserInitialization(pDetect);

  // Initialize the default Hadrontherapy geometry
  geometryController->SetGeometry("default");

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

  // Interaction data
  new HadrontherapyInteractionParameters(pDetect);

#ifdef G4VIS_USE
  // Visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager -> Initialize();
#endif 

G4UImanager* UI = G4UImanager::GetUIpointer();      
  
 if (argc!=1)   // batch mode
   {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);    
   }
 
 else  // interactive mode : define visualization UI terminal
   {
     G4UIsession* session = 0;
     
     // If the enviroment variable for the TCSH terminal is active, it is used and the
     // defaultMacro.mac file is executed
#if defined(G4UI_USE_TCSH)
     session = new G4UIterminal(new G4UItcsh);      
     UI->ApplyCommand("/control/execute defaultMacro.mac");  

     // Alternatively (if G4UI_USE_TCSH is not defined)  the program search for the
     // G$UI_USE_QT variable. It starts a graphical user interface based on the QT libraries
     // In this case a gui.mac file is executed
#elif defined(G4UI_USE_QT)
     session = new G4UIQt(argc,argv);
     UI->ApplyCommand("/control/execute gui.mac");      
     
     // As final option, the simpler user interface terminal is opened
#else
     session = new G4UIterminal();
#endif
      
      //#ifdef G4VIS_USE
      //UI->ApplyCommand("/control/execute vis.mac");     
      //#endif
      
	session->SessionStart();
	delete session;
   }
 matrix -> TotalEnergyDeposit();
 
#ifdef ANALYSIS_USE
  analysis -> finish();
#endif

  // Job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  
  delete geometryMessenger;
  delete geometryController;
  delete runManager;
  return 0;
}
