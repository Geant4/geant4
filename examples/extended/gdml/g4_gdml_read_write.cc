//
// -*- C++ -*-
// $Id: g4_gdml_read_write.cc,v 1.2 2004-12-06 18:13:04 gcosmo Exp $
//
// Author: Radovan Chytracek 2000 - 2004
//         Witek Pokorski
//
#include <stdexcept>

#include "globals.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4TransportationManager.hh"

#ifdef G4VIS_USE
#include "g4rwVisManager.hh"
#endif

#include "g4rwDetectorConstruction.hh"
#include "g4rwPhysicsList.hh"
#include "g4rwPrimaryGeneratorAction.hh"

//g4 writer
#include "WriterG4/G4GDMLWriter.h"


int main()
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;
  //
#ifdef G4VIS_USE
  G4VisManager* visManager = new gogdmlVisManager;
  visManager->Initialize();
#endif

  // set mandatory initialization classes
  runManager->SetUserInitialization(new gogdmlDetectorConstruction);
  runManager->SetUserInitialization(new gogdmlPhysicsList);

  // set mandatory user action class
  runManager->SetUserAction(new gogdmlPrimaryGeneratorAction);

  // initialize G4 kernel
  runManager->Initialize();  

  // scanning geometry tree
  G4VPhysicalVolume* g4wv =
    G4TransportationManager::GetTransportationManager()->
    GetNavigatorForTracking()->GetWorldVolume();
  
  G4GDMLWriter g4writer("schema/gdml_2.0.xsd", "g4writertest.gdml");
  
  try
  {
    g4writer.DumpGeometryInfo(g4wv);
  }
  catch(std::logic_error &lerr)
  {
    std::cout << "Caught an exception: " 
              << lerr.what () << std::endl;
  }

  
  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();

  // G4UIterminal is a (dumb) terminal.
  G4UIsession * session = new G4UIterminal(new G4UItcsh);
  
  UI->ApplyCommand("/control/execute vis.mac"); 

  UI->ApplyCommand("/run/verbose 0");
  UI->ApplyCommand("/event/verbose 0");
  UI->ApplyCommand("/tracking/verbose 1");

  // start a run
  int numberOfEvent = 10;
  runManager->BeamOn(numberOfEvent);

  session->SessionStart();

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}


