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
// $Id: thyroid.cc,v 1.1 2002-06-12 14:03:22 francy Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------
//                 GEANT 4 - Thyroid example
// --------------------------------------------------------------
//
// Code developed by:
// S. Agostinelli, F. Foppiano, S. Garelli and M. Tropeano
//
// Brachytherapy simulates the dose deposition in a cubic (30*cm)
// water phantom for a Ir-192 MicroSelectron High Dose Rate
// brachytherapy source.
//
// Simplified gamma generation is used.
// Source axis is oriented along Z axis. 
// Voxel data on the X-Z plane is output to ASCII file
// "Brachy.out".
//
// For information related to this code contact the developers.
//
// --------------------------------------------------------------  

#include "ThyroidEventAction.hh"
#include "ThyroidDetectorConstruction.hh"
#include "ThyroidPhysicsList.hh"
#include "ThyroidPrimaryGeneratorAction.hh"
#include "ThyroidSD.hh"
//#include "ThyroidSteppingVerbose.hh"

#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4VIS_USE
#include "ThyroidVisManager.hh"
#endif

//#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
int main(int argc ,char ** argv)
{
 // Number of generated photons
  //G4int NumberOfEvents = 10;

  //my Verbose output class
//  G4VSteppingVerbose::SetInstance(new ThyroidSteppingVerbose);


 // Construct the default run manager
 G4RunManager* pRunManager = new G4RunManager;

#ifdef G4VIS_USE
  // visualization manager
 G4VisManager* visManager = new ThyroidVisManager;
 visManager->Initialize();
#endif


 // Set mandatory initialization classes
 ThyroidDetectorConstruction *pDetectorConstruction;
 pRunManager->SetUserInitialization(pDetectorConstruction = new
ThyroidDetectorConstruction); 
pRunManager->SetUserInitialization(new ThyroidPhysicsList);

 // Set mandatory user action class
 pRunManager->SetUserAction(new ThyroidPrimaryGeneratorAction);
 pRunManager->SetUserAction(new ThyroidEventAction);


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
 



 // Initialize G4 kernel
 pRunManager->Initialize();


 // Get the pointer to the UI manager and set verbosities
 G4UImanager* UI = G4UImanager::GetUIpointer();
 UI->ApplyCommand("/run/verbose 0");
 UI->ApplyCommand("/event/verbose 0");
 UI->ApplyCommand("/tracking/verbose 1");

 // pRunManager->BeamOn(NumberOfEvents);

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
 
  
	



	G4std::ofstream ofs;

	// Output voxel data to text file
	// Format = x coord [mm] <tab> y coord [mm]<tab> z coord [mm] <tab> edep [MeV] <eol>
	/*
            ofs.open("Thyroid.out");
		{ 
		  G4double Value_X = pDetectorConstruction-> m_theDx;
                  G4double Value_Y= pDetectorConstruction-> m_theDy;
                  G4double Value_Z = pDetectorConstruction-> m_theDz;
                  ofs << Value_X << '\t' << Value_Y << '\t' << Value_Z << G4endl;
			
		ofs.close();

 }
	*/

#ifdef G4VIS_USE
  delete visManager;
#endif
  


 // Job termination
 delete pRunManager;

 return 0;
}












