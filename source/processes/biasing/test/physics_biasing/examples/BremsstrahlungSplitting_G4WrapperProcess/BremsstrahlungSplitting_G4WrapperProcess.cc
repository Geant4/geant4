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
// $Id: BremsstrahlungSplitting_G4WrapperProcess.cc,v 1.1 2007-05-24 21:57:02 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation.
//

#include "ConfigData.hh"
#include "ConfigDataMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UItcsh.hh"
#include "G4UIterminal.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4UI_USE_WIN32
#include "G4UIWin32.hh"
#endif

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

int main(int argc,char** argv) 
{
  CLHEP::HepRandom::setTheEngine(new CLHEP::MTwistEngine);
     
  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // Mandatory initialization classes
  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new PhysicsList);

  // Configuration data messenger
  ConfigDataMessenger* configDataMsgr = new ConfigDataMessenger();

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
    
  // User action classes
  runManager->SetUserAction(new PrimaryGeneratorAction());
  runManager->SetUserAction(new RunAction);
  
  G4UIsession* session=0;
  G4UImanager* UI = G4UImanager::GetUIpointer();  
  
  if (argc==1) {
    // Define UI session for interactive mode.
    // G4UIterminal is a (dumb) terminal.
    
#if defined(G4UI_USE_XM)
    session = new G4UIXm(argc,argv);
#elif defined(G4UI_USE_WIN32)
    session = new G4UIWin32();
#elif defined(G4UI_USE_TCSH)
    session = new G4UIterminal(new G4UItcsh);      
#else
    session = new G4UIterminal();
#endif
    session->SessionStart();
    delete session;
  }
  else { 
    // Batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    G4String job = argv[2];
    G4String geometryFile = argv[3];
    G4String scoringFile = argv[4];
    G4String physicsFile = argv[5];
    G4String cutsFile = argv[6];
    G4String biasingFile = argv[7];
    
    G4cout<<"Executing "<<command+fileName<<" "<<", for geometry "<<geometryFile<<", for scoring "<<scoringFile<<", for physics "<<physicsFile<<", for cuts "<<cutsFile<<", and biasing "<<biasingFile<<G4endl;
    
    
    G4String outputDir = "Results_"+job;
    G4String outputId = geometryFile+":"+scoringFile+":"+physicsFile+":"+cutsFile+":"+biasingFile;
    
    ConfigData::SetOutputDirectory(outputDir);
    ConfigData::SetOutputId(outputId);
    
    // UI->ApplyCommand("/control/shell rm -r " + outputDir);
    UI->ApplyCommand("/control/shell mkdir -p " + outputDir);
    
    UI->ApplyCommand("/control/alias geometry Macros/Geometry/" + geometryFile+".mac");
    UI->ApplyCommand("/control/alias scoring Macros/Scoring/" + scoringFile+".mac");
    UI->ApplyCommand("/control/alias physics Macros/Physics/" + physicsFile+".mac");
    UI->ApplyCommand("/control/alias cuts Macros/Cuts/" + cutsFile+".mac");
    UI->ApplyCommand("/control/alias biasing Macros/Biasing/" + biasingFile+".mac");
    UI->ApplyCommand(command+fileName);
  }
  
  delete visManager;
  delete runManager;
  delete configDataMsgr;

  return 0;
}
