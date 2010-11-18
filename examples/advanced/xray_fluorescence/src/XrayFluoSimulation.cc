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
// $Id: XrayFluoSimulation.cc
//
// Author: Elena Guardincerri 
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
// 24 Ago 2002 Splitted in a separet class  Alfonso Mantero
//
// -------------------------------------------------------------------
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif
#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif
#include "XrayFluoDetectorConstruction.hh"
#include "XrayFluoPlaneDetectorConstruction.hh"
#include "XrayFluoMercuryDetectorConstruction.hh"
#include "XrayFluoPhysicsList.hh"
#include "XrayFluoPrimaryGeneratorAction.hh"
#include "XrayFluoPlanePrimaryGeneratorAction.hh"
#include "XrayFluoMercuryPrimaryGeneratorAction.hh"
#include "XrayFluoRunAction.hh"
#include "XrayFluoEventAction.hh"
#include "XrayFluoSteppingAction.hh"
#include "XrayFluoSteppingVerbose.hh"
#include "XrayFluoSimulation.hh"
#ifdef G4ANALYSIS_USE
#include "XrayFluoAnalysisManager.hh"
#endif
using namespace CLHEP;

XrayFluoSimulation::XrayFluoSimulation(G4int seed):dir(seed)

{ }

XrayFluoSimulation::~XrayFluoSimulation()

{ }

void XrayFluoSimulation::RunSimulation(int argc,char* argv[])
{

  // choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(dir);

  //XrayFluo Verbose output class
  G4VSteppingVerbose::SetInstance(new XrayFluoSteppingVerbose);
     
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization 

  XrayFluoPhysicsList* xrayList = 0;

  // chosing Geometry setup

  G4int geometryNumber; 

  if (argc == 3){
    geometryNumber =  atoi(argv[2]);

  }
  while ( (geometryNumber != 1) && (geometryNumber !=2) && (geometryNumber !=3) && (geometryNumber !=4)) {
    G4cout << "Please Select Simulation Geometrical Set-Up: "<< G4endl;
    G4cout << "1 - Test Beam" << G4endl;
    G4cout << "2 - Infinite Plane" << G4endl;
    G4cout << "3 - Planet and Sun"<< G4endl;
    G4cout << "4 - Phase-Space Production"<< G4endl;


    G4cin >> geometryNumber;
  }

  XrayFluoDetectorConstruction* testBeamDetector = 0;
  XrayFluoPlaneDetectorConstruction* planeDetector = 0;
  XrayFluoMercuryDetectorConstruction* mercuryDetector = 0;



  if (geometryNumber == 1 || geometryNumber == 4) {
    testBeamDetector = XrayFluoDetectorConstruction::GetInstance();
    if (geometryNumber == 4) {
      testBeamDetector->PhaseSpaceOn();
    }
    runManager->SetUserInitialization(testBeamDetector);
    xrayList = new   XrayFluoPhysicsList(testBeamDetector);
  }
  else if (geometryNumber == 2) {
    planeDetector = XrayFluoPlaneDetectorConstruction::GetInstance();
    runManager->SetUserInitialization(planeDetector);
    xrayList = new   XrayFluoPhysicsList(planeDetector);
  }
  else if (geometryNumber == 3) {
    mercuryDetector = XrayFluoMercuryDetectorConstruction::GetInstance();
    runManager->SetUserInitialization(mercuryDetector);
    xrayList = new   XrayFluoPhysicsList(mercuryDetector);
  }


  runManager->SetUserInitialization(xrayList);
  
#ifdef G4VIS_USE
  //visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

#ifdef G4ANALYSIS_USE
  // set analysis to have the messenger running...
  XrayFluoAnalysisManager* analysis = XrayFluoAnalysisManager::getInstance();
#endif
  XrayFluoEventAction* eventAction = 0;
  XrayFluoRunAction* runAction = new XrayFluoRunAction();
  XrayFluoSteppingAction* stepAction = new XrayFluoSteppingAction();


  //Selecting the PrimaryGenerator depending upon Geometry setup selected

  if (geometryNumber == 1 || geometryNumber == 4) {
    if (geometryNumber == 4) {
#ifdef G4ANALYSIS_USE
     analysis->PhaseSpaceOn();
     analysis->CreatePersistency();
#endif
    }
    eventAction = new XrayFluoEventAction(testBeamDetector);
    runManager->SetUserAction(new XrayFluoPrimaryGeneratorAction(testBeamDetector));
  }

  else if (geometryNumber == 2) {
    eventAction = new XrayFluoEventAction(planeDetector);
    runManager->SetUserAction(new XrayFluoPlanePrimaryGeneratorAction(planeDetector));
  }

  else if (geometryNumber == 3) {
   stepAction->SetMercuryFlag(true);
   eventAction = new XrayFluoEventAction(mercuryDetector);
   runManager->SetUserAction(new XrayFluoMercuryPrimaryGeneratorAction(mercuryDetector));
      }
 
  runManager->SetUserAction(eventAction); 
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(stepAction);
  
  //Initialize G4 kernel
  runManager->Initialize();
  
  // get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  if (argc == 1)   // Define UI session for interactive mode.
    {
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute initInter.mac");     
#endif
      if (ui->IsGUI())
        UImanager->ApplyCommand("/control/execute gui.mac");     
      ui->SessionStart();
      delete ui;
#endif
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  
  // job termination
  
#ifdef G4VIS_USE
  delete visManager;
  G4cout << "visManager deleted"<< G4endl;
#endif

//   if (testBeamDetector) delete testBeamDetector;
//   if (planeDetector) delete planeDetector;
//   if (mercuryDetector) delete mercuryDetector;

  
  delete runManager;
  
   

}
