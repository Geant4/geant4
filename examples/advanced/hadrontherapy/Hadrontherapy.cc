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
// $Id: Hadrontherapy.cc
//
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
//
// ****************  Hadrontherapy  ****************************************
// Hadrontherapy simulates a general transport beam line dedicated to the 
// irradiation of toumors with hadron beams.
// All the elements of a typical hadron beam line (collimator, 
// scattering system, range shifter, etc.) are simulated.
// Positions, dimensions and materials of such element can be changed by the users.
// Actually only proton beams can be simulated.
// All the characteristics of the incident beam can be changed.
// Two typical detectors commonly used in the hadrontherapy
// field are simulated: the Markus ionization chamber for the 
// reconstruction of the depth dose distributions,
// and a gafchromic film for the reconstruction of the 
// lateral dose distributions.
// **************************************************************************

#include <fstream>
#include <iomanip>
#include <iostream>
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "HadrontherapyVisManager.hh"
#endif

#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapyRunAction.hh"
#include "HadrontherapyEventAction.hh"
#include "HadrontherapySteppingAction.hh"
// -----------------------------------------------------------------------
int main(int argc,char** argv) {

  //***************************
  // choose the Random engine
  //***************************

  HepRandom::setTheEngine(new RanecuEngine);
  G4int seed = time(NULL);
  HepRandom::setTheSeed(seed);
    
  //***********************************
  // Construct the default run manager
  //***********************************
  
  G4RunManager * runManager = new G4RunManager;
  
  //***************************************
  // set mandatory initialization classes
  //***************************************

  HadrontherapyDetectorConstruction* detector;
  detector = new HadrontherapyDetectorConstruction;
  runManager -> SetUserInitialization(detector);
  runManager -> SetUserInitialization(new HadrontherapyPhysicsList(detector));

  //***********************************************
  // Set the visualization if you chose to have it!
  //***********************************************

#ifdef G4VIS_USE
  G4VisManager* visManager = new HadrontherapyVisManager;
  visManager -> Initialize();
#endif 

  //**********************************
  // set mandatory user action class
  //********************************

  runManager -> SetUserAction(new HadrontherapyPrimaryGeneratorAction( detector ));
 
  //****************************************
  // set the optional user action classes
  //***************************************

  HadrontherapyRunAction* runaction = new HadrontherapyRunAction;
  runManager -> SetUserAction(runaction);
 
  HadrontherapyEventAction* eventaction = new HadrontherapyEventAction( runaction );
  runManager -> SetUserAction(eventaction);
  
  HadrontherapySteppingAction* steppingaction = new HadrontherapySteppingAction( eventaction );
  runManager -> SetUserAction(steppingaction);    
 
  //*********************
  // Initialize G4 kernel
  //*********************

  runManager -> Initialize();

  //***********************************************
  // get the pointer to the User Interface manager 
  //***********************************************

  G4UImanager* UI = G4UImanager::GetUIpointer();

  //*******************************************************************
  //Define  UI terminal for interactive mode (wait command from keyboard
  //or for batch mode but reading a macro file
  //********************************************************************

  G4UIsession* session = 0;
  
  if (argc==1)   // Define UI session for interactive mode.
    {
                        
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


  if (session)   // Define UI session for interactive mode.
    {
      
      UI->ApplyCommand("/control/execute defaultMacro.mac");    
      session -> SessionStart();
      delete session;
    }

  else           // Batch mode

    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }
 
  //******************* 
  // job termination
  //*******************

#ifdef G4VIS_USE
  delete visManager;
#endif  

  delete runManager;

  return 0;
}

