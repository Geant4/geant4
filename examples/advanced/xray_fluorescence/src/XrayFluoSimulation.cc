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
// 18 Jan 2011 adapted to new deexcitation design
//
// -------------------------------------------------------------------

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

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
#include "XrayFluoSimulation.hh"

#include "XrayFluoActionInitializer.hh"
#include "XrayFluoAnalysisManager.hh"


XrayFluoSimulation::XrayFluoSimulation(G4int seed):dir(seed)
{;}


XrayFluoSimulation::~XrayFluoSimulation()
{;}

void XrayFluoSimulation::RunSimulation(int argc,char* argv[])
{
  // choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4Random::setTheSeed(dir);
   
#ifdef G4MULTITHREADED
  // Construct the default run manager
  G4MTRunManager* runManager = new G4MTRunManager();
  G4cout << "Using the MT Run Manager (G4MULTITHREADED=ON)" << G4endl; 
#else
  G4RunManager * runManager = new G4RunManager;
  G4cout << "Using the sequential Run Manager" << G4endl;
#endif


  // chosing Geometry setup
  G4int geometryNumber = 0; 

  if (argc == 3){
    geometryNumber =  atoi(argv[2]);

  }
  while ( (geometryNumber < 1) || (geometryNumber >4)) {
    G4cout << "Please Select Simulation Geometrical Set-Up: "<< G4endl;
    G4cout << "1 - Test Beam" << G4endl;
    G4cout << "2 - Infinite Plane" << G4endl;
    G4cout << "3 - Planet and Sun"<< G4endl;
    G4cout << "4 - Phase-Space Production"<< G4endl;

    G4cin >> geometryNumber;
  }

  // set analysis to have the messenger running...
  XrayFluoAnalysisManager* analysis = 
    XrayFluoAnalysisManager::getInstance();

  // set mandatory initialization 

  //Initialize geometry
  if (geometryNumber == 1 || geometryNumber == 4) {
    XrayFluoDetectorConstruction* testBeamDetector 
      = XrayFluoDetectorConstruction::GetInstance();
    if (geometryNumber == 4) {
      testBeamDetector->PhaseSpaceOn();
      analysis->PhaseSpaceOn();
    }
    runManager->SetUserInitialization(testBeamDetector);
  }
  else if (geometryNumber == 2) {
    XrayFluoPlaneDetectorConstruction* 
      planeDetector = XrayFluoPlaneDetectorConstruction::GetInstance();
    runManager->SetUserInitialization(planeDetector);
  }
  else if (geometryNumber == 3) {
    XrayFluoMercuryDetectorConstruction* mercuryDetector = 
      XrayFluoMercuryDetectorConstruction::GetInstance();
    runManager->SetUserInitialization(mercuryDetector);
  }


  //Initialize physics
  runManager->SetUserInitialization(new   XrayFluoPhysicsList()); 

  //Initialize actions
  runManager->SetUserInitialization
    (new XrayFluoActionInitializer(geometryNumber));

#ifdef G4VIS_USE
  //visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
  
  // get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  if (getenv("G4VIS_USE")) {
    UImanager->ApplyCommand("/control/execute vis.mac");
  }     

  if (argc == 1)   // Define UI session for interactive mode.
    {
      UImanager->ApplyCommand("/control/execute initInter.mac");     
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);        
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

  
  delete runManager;
  
  
}
