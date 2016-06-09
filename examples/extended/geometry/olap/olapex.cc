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
// $Id: olapex.cc,v 1.3 2006/06/29 17:21:55 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// 
// --------------------------------------------------------------
// GEANT 4 - OLAP, Debugging Tool for Ovelapping Geometries
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
// the detector-construction
#include "RandomDetector.hh"

// Geant4
#include "G4RunManager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

// this module
#include "OlapGenerator.hh"
#include "OlapDetConstr.hh"
#include "OlapPhysicsList.hh"
#include "OlapManager.hh"
#include "OlapRunAction.hh"
#include "OlapEventAction.hh"
#include "OlapSteppingAction.hh"
#include "G4GeoNav.hh"

#ifdef G4VIS_USE
  #include "G4VisExecutive.hh"
#endif

#ifdef debug
 #include "G4LogicalVolumeStore.hh"
 #include "G4PhysicalVolumeStore.hh"
#endif


int main(int argc,char** argv) {


  G4RunManager * runManager = new G4RunManager;
  
  // particle-generator must be the first useraction!!!
  runManager->SetUserAction(new OlapGenerator);
  
  // put the user-geometry here
  G4VUserDetectorConstruction * origGeom = new RandomDetector(0);
  OlapDetConstr * olapGeom = new OlapDetConstr(origGeom);
  runManager->SetUserInitialization(olapGeom);
  runManager->SetUserInitialization(new OlapPhysicsList); 
  
  // initialize G4
  runManager->Initialize();
  
  // initialize overlap-detection
  OlapManager * olap = OlapManager::GetOlapManager();
  olap->GotoLV(olap->GetFullWorld()->GetLogicalVolume()->GetName());
  
#ifdef G4VIS_USE
  //Instantiating VisManager
  G4VisManager *visManager = new G4VisExecutive;
  visManager -> Initialize();
#endif

  //User action classes.
  OlapRunAction * olapRunAction = new OlapRunAction;
  OlapEventAction * olapEventAction = new OlapEventAction(olapRunAction);
  OlapSteppingAction * olapSteppingAction = new OlapSteppingAction(olapEventAction);
  runManager->SetUserAction(olapRunAction);
  runManager->SetUserAction(olapEventAction);
  runManager->SetUserAction(olapSteppingAction);

  

#ifdef debug
  //Get number of logical and physical volumes
  G4LogicalVolumeStore*  logstore = G4LogicalVolumeStore::GetInstance();
  G4PhysicalVolumeStore* phystore = G4PhysicalVolumeStore::GetInstance();
  G4cout << endl
       << "Number of Logical volumes     " << logstore->entries() << endl
       << "====================================" << endl
       << "Number of Physical volumes    " << phystore->entries() << endl
       << "====================================" << endl << endl;
#endif
   
   

  //Get the pointer for the User Interface maager
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if (argc == 1) {
    G4VUIshell * s = new G4UItcsh("Idle> ",10);
    G4UIsession* session = new G4UIterminal(s);
    
    G4String command = "/control/execute gui.mac";
    UI->ApplyCommand(command);        
    session->SessionStart();

    delete session;
  }
  //Batch mode
  else {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);    
  }


#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}
