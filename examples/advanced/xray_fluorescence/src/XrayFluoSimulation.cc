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
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif
#include "Randomize.hh"
#ifdef G4VIS_USE
#include "XrayFluoVisManager.hh"
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

XrayFluoSimulation::XrayFluoSimulation(G4int seed):dir(seed)

{ }

XrayFluoSimulation::~XrayFluoSimulation()

{ }

void XrayFluoSimulation::RunSimulation(int argc,char* argv[])
{

  // choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  HepRandom::setTheSeed(dir);

  //XrayFluo Verbose output class
  G4VSteppingVerbose::SetInstance(new XrayFluoSteppingVerbose);
     
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization 

  XrayFluoPhysicsList* xrayList = 0;

  // chosing Geometry setup

  G4int GeometryNumber;


  if (argc == 3){
    GeometryNumber =  atoi(argv[2]);
  }
  while ( (GeometryNumber != 1) && (GeometryNumber !=2) && (GeometryNumber !=3) ) {
    G4cout << "Please Select Simulation Geometrical Set-Up: "<< G4endl;
    G4cout << "1 - Test Beam" << G4endl;
    G4cout << "2 - Infinite Plane" << G4endl;
    G4cout << "3 - Planet and Sun (beta)"<< G4endl;
    G4cin >> GeometryNumber;
  }

  XrayFluoDetectorConstruction* testBeamDetector = 0;
  XrayFluoPlaneDetectorConstruction* planeDetector = 0;
  XrayFluoMercuryDetectorConstruction* mercuryDetector = 0;

  if (GeometryNumber == 1) {
    testBeamDetector = XrayFluoDetectorConstruction::GetInstance();
    runManager->SetUserInitialization(testBeamDetector);
    xrayList = new   XrayFluoPhysicsList(testBeamDetector);
  }
  else if (GeometryNumber == 2) {
    planeDetector = XrayFluoPlaneDetectorConstruction::GetInstance();
    runManager->SetUserInitialization(planeDetector);
    xrayList = new   XrayFluoPhysicsList(planeDetector);
  }
  else if (GeometryNumber == 3) {
    mercuryDetector = XrayFluoMercuryDetectorConstruction::GetInstance();
    runManager->SetUserInitialization(mercuryDetector);
    xrayList = new   XrayFluoPhysicsList(mercuryDetector);
  }

  runManager->SetUserInitialization(xrayList);
  
  G4UIsession* session=0;
  
  if (argc==1)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
#ifdef G4UI_USE_XM
      session = new G4UIXm(argc,argv);
#else           
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);
#else
      session = new G4UIterminal();
#endif
#endif
    }


  
#ifdef G4VIS_USE
  //visualization manager
  G4VisManager* visManager = new XrayFluoVisManager;
  visManager->Initialize();
#endif

  // set analysis to have the messenger running...
  // XrayFluoAnalysisManager* analysis =

#ifdef G4ANALYSIS_USE
  XrayFluoAnalysisManager::getInstance();
#endif

  XrayFluoEventAction* eventAction = 0;
  XrayFluoRunAction* runAction = new XrayFluoRunAction();
  XrayFluoSteppingAction* stepAction = new XrayFluoSteppingAction();


  //Selecting the PrimaryGenerator depending upon Geometry setup selected

  if (GeometryNumber == 1) {
    eventAction = new XrayFluoEventAction(testBeamDetector);
    runManager->SetUserAction(new XrayFluoPrimaryGeneratorAction(testBeamDetector));

  }
  else if (GeometryNumber == 2) {
    eventAction = new XrayFluoEventAction(planeDetector);
    runManager->SetUserAction(new XrayFluoPlanePrimaryGeneratorAction(planeDetector));
  }
  else if (GeometryNumber == 3) {
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
  G4UImanager* UI = G4UImanager::GetUIpointer();
  
  if (session)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
      UI->ApplyCommand("/control/execute initInter.mac");    
#ifdef G4UI_USE_XM
      // Customize the G4UIXm menubar with a macro file :
      UI->ApplyCommand("/control/execute gui.mac");
#endif
      session->SessionStart();
      delete session;
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }
  
  // job termination
  
#ifdef G4VIS_USE
  delete visManager;
  G4cout << "visManager deleted"<< G4endl;
#endif
  
  delete runManager;
  
  if (testBeamDetector) delete testBeamDetector;
  if (planeDetector) delete planeDetector;
  if (mercuryDetector) delete mercuryDetector;
   

}
