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
// --------------------------------------------------------------
//                 GEANT 4 - Brachytherapy example
// --------------------------------------------------------------
//
// Code developed by:
// S. Agostinelli, F. Foppiano, S. Garelli , M. Tropeano,S.Guatelli
//
// Brachytherapy simulates the dose deposition in a cubic (30*cm)
// water phantom for a IsoSeed I-125 Bebig
// brachytherapy source.
//
// Simplified gamma generation is used.
// Source axis is oriented along Z axis. 
// Voxel data on the X-Z plan
//e is output to 
// "Brachy3.hbk".


// $Id: Brachy.cc,v 1.12 2002-06-25 13:12:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4VIS_USE
#include "BrachyVisManager.hh"
#endif

#include "BrachyEventAction.hh"
#include "BrachyDetectorConstruction.hh"
#include "BrachyPhysicsList.hh"
#include "BrachyPrimaryGeneratorAction.hh"
#include "BrachyPhantomSD.hh"

#include "G4SDManager.hh"

#include"BrachyRunAction.hh"

#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"



//Interactive mode//

int main(int argc ,char ** argv)

{ 
  // semi diversi li ha usati per l'errore
  

  //HepRandom::setTheSeed(16520);

  HepRandom::setTheEngine(new RanecuEngine);
 
  
  // Construct the default run manager
  G4RunManager* pRunManager = new G4RunManager;


  // Set mandatory initialization classes
  G4String SDName = "Phantom";
 
  BrachyDetectorConstruction *pDetectorConstruction=new BrachyDetectorConstruction(SDName);
  pRunManager->SetUserInitialization(pDetectorConstruction) ;
 
  pRunManager->SetUserInitialization(new BrachyPhysicsList);

  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new BrachyVisManager;
  visManager->Initialize();
#endif

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
 


  BrachyEventAction *pEventAction=new BrachyEventAction(SDName);
 
  pRunManager->SetUserAction(pEventAction );
  pRunManager->SetUserAction(new BrachyPrimaryGeneratorAction);
 
  BrachyRunAction *pRunAction=new BrachyRunAction();
  pRunManager->SetUserAction(pRunAction);
 
 
  
//Initialize G4 kernel
  pRunManager->Initialize();

  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  
  UI->ApplyCommand("/run/verbose 0");
  UI->ApplyCommand("/event/verbose 0");
  UI->ApplyCommand("/tracking/verbose 0");

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
  
  //int numberOfEvent = 10000;
  //pRunManager->BeamOn(numberOfEvent);
  				
 	   
	 
  // Job termination


  
#ifdef G4VIS_USE
  delete visManager;
#endif
  
   

  delete pRunManager;

  return 0;
}








