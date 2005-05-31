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
// $Id: PhotIn.cc,v 1.2 2005-05-31 15:23:01 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - example/advanced/Photon_inefficiency 
//
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------

#define debug

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "PhotInVisManager.hh"
#endif
#include "PhotInDetectorConstruction.hh"
#include "PhotInDetectorMessenger.hh"
#include "PhotInPhysicsList.hh"
#include "PhotInPrimaryGeneratorAction.hh"
#include "PhotInRunAction.hh"
#include "PhotInEventAction.hh"
#include "PhotInStackingAction.hh"

int main(int argc,char** argv)
{
#ifdef debug
  G4cout<<"PhotIn: Example is started"<<G4endl;
#endif

 // Construct the default run manager
 G4RunManager* runManager = new G4RunManager; // Del? M.K.

#ifdef debug
  G4cout<<"PhotIn: Run manager is created"<<G4endl;
#endif

 // set mandatory initialization classes
		G4double hx=.5*m;
		G4double hy=.5*m;
		G4double hz= m;
  PhotInDetectorConstruction* detector = new PhotInDetectorConstruction(hx,hy,hz); // Del?
  //detector->SetSectionHalfDimensions(0.5*m,0.5*m,m);// This is how it could be changed(-)

#ifdef debug
  G4cout<<"PhotIn: Detector is created with hx="<<hx<<", hy="<<hy<<", hz="<<hz<<G4endl;
#endif

		//PhotInDetectorMessenger* detectorMessenger=new PhotInDetectorMessenger(detector);//Del?
  new PhotInDetectorMessenger(detector);// Del?

#ifdef debug
  G4cout<<"PhotIn: Detector Messenger is constructed"<<G4endl;
#endif

  runManager->SetUserInitialization(detector);

#ifdef debug
  G4cout<<"PhotIn: Detector Messenger is constructed"<<G4endl;
#endif

  PhotInPhysicsList* physList = new PhotInPhysicsList; // Del? M.K.

#ifdef debug
  G4cout<<"PhotIn: Physics List is constructed"<<G4endl;
#endif

  runManager->SetUserInitialization(physList);

#ifdef debug
  G4cout<<"PhotIn: Physics List is transfered to Run Manager"<<G4endl;
#endif
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new PhotInVisManager;
  visManager->Initialize();
#ifdef debug
  G4cout<<"PhotIn: Visualization Manager is constructed and initialized"<<G4endl;
#endif

#endif
    
  // set user action classes
  PhotInPrimaryGeneratorAction* primaryGeneration = new PhotInPrimaryGeneratorAction;
  runManager->SetUserAction(primaryGeneration);

#ifdef debug
  G4cout<<"PhotIn: Primary Generator is constructed"<<G4endl;
#endif

  //PhotInPrimaryGeneratorAction* gen = (PhotInPrimaryGeneratorAction*)
  //                       G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  //if(gen) gen->SetDetector(detector);
  //else G4cout<<"---Warning--- PhotIn: PhotInPrimaryGeneratorAction* ="<<gen<<G4endl;
  // In the example detectorMessenger is created only for DetectorConstruction to change it

  primaryGeneration->SetDetector(detector);

#ifdef debug
  G4cout<<"PhotIn: Detector is transfered to the primary generator"<<G4endl;
#endif

  runManager->SetUserAction(new PhotInRunAction);

#ifdef debug
  G4cout<<"PhotIn: Run is constructed"<<G4endl;
#endif

  runManager->SetUserAction(new PhotInEventAction);

#ifdef debug
  G4cout<<"PhotIn: Event List is constructed"<<G4endl;
#endif

  runManager->SetUserAction(new PhotInStackingAction);
  
#ifdef debug
  G4cout<<"PhotIn: Stacking Actions are constructed"<<G4endl;
#endif

  // Initialize G4 kernel
  // runManager->Initialize(); // (?) (M.K.)
     
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();

#ifdef debug
  G4cout<<"PhotIn: Start UI initialization"<<G4endl;
#endif  

  if (argc==1)   // Define UI session for interactive mode.
  {
    G4UIsession* session=0;
  
    // G4UIterminal is a (dumb) terminal.
#ifdef G4UI_USE_TCSH
    session = new G4UIterminal(new G4UItcsh);      
#else
    session = new G4UIterminal();
#endif    

#ifdef debug
    G4cout<<"PhotIn: UITerminal is constructed"<<G4endl;
#endif
      
    UI->ApplyCommand("/control/execute vis.mac");    

#ifdef debug
    G4cout<<"PhotIn: First UI comand /control/execute vis.mac is sent"<<G4endl;
#endif
    session->SessionStart();

#ifdef debug
    G4cout<<"PhotIn: UI session is started"<<G4endl;
#endif

    delete session;

#ifdef debug
    G4cout<<"PhotIn: UI session is finished"<<G4endl;
#endif

  }
  else           // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];

#ifdef debug
    G4cout<<"PhotIn: Batch mode /control/execute "<<argv[1]<<" is started"<<G4endl;
#endif

    UI->ApplyCommand(command+fileName);
  }

  // job termination
#ifdef G4VIS_USE
  delete visManager;

#ifdef debug
  G4cout<<"PhotIn: Visualization Manager is deleted"<<G4endl;
#endif

#endif

  //delete detectorMessenger; // ?

#ifdef debug
  G4cout<<"PhotIn: Detector Messenger is deleted"<<G4endl;
#endif

  delete runManager;

#ifdef debug
  G4cout<<"PhotIn: Run Manager is deleted"<<G4endl;
#endif

  return 0;
}

