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
// $Id: olapex.cc,v 1.1 2002-06-04 07:40:18 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  #include "OlapVisManager.hh"
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
  //Instantiating OlapVisManager
  OlapVisManager *visManager = new OlapVisManager;
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
