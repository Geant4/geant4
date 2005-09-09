//
// File name:     RadmonApplication.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplication.cc,v 1.1 2005-09-09 08:26:54 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonApplication.hh"
#include "RadmonApplicationOptions.hh"

#include "RadmonDetectorLayout.hh"
#include "RadmonDetectorConstruction.hh"
#include "RadmonDetectorLabelledEntitiesConstructorsFactory.hh"
#include "RadmonDetectorMessenger.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#ifdef    G4VIS_USE
 #include "G4VisExecutive.hh"
#endif /* G4VIS_USE */



                                                RadmonApplication :: RadmonApplication(const RadmonApplicationOptions & options)
:
 valid(false),
 runManager(0),
 layout(0),
 factory(0),
 #ifdef    G4VIS_USE
  visManager(0),
 #endif /* G4VIS_USE */
 uiManager(0),
 messenger(0),
 session(0)
{
 // Construct the default run manager
 runManager=new G4RunManager();
 
 if (runManager==0)
 {
  G4cerr << options.ApplicationName() << ": Run manager not allocated." << G4endl;
  return;
 }
 
 
 // Construct the detector layout
 layout=new RadmonDetectorLayout;
 
 if (layout==0)
 {
  G4cerr << options.ApplicationName() << ": Layout not allocated." << G4endl;
  return;
 }
 
 
 // Construct the detectors factory
 factory=new RadmonDetectorLabelledEntitiesConstructorsFactory;
 
 if (factory==0)
 {
  G4cerr << options.ApplicationName() << ": Factory not allocated." << G4endl;
  return;
 }
 
 
 // Construct the entity constructors
 if (!CreateEntityConstructors())
 {
  G4cerr << options.ApplicationName() << ": Entity constructors not allocated." << G4endl;
  return;
 }
 
 
 // Construct the detector construction
 RadmonDetectorConstruction * detectorConstruction(new RadmonDetectorConstruction(layout, factory));
 
 if (detectorConstruction==0)
 {
  G4cerr << options.ApplicationName() << ": Detector construction not allocated." << G4endl;
  return;
 }
 
 
 // The factory will is owned by the detectorConstruction
 factory=0; 
 
 runManager->SetUserInitialization(detectorConstruction);
 
 // Set mandatory initialization classes
 // runManager->SetUserInitialization(new ExN01PhysicsList);
         
 // Set mandatory user action class
 // runManager->SetUserAction(new ExN01PrimaryGeneratorAction);
 
            
 // Initialize the run manager
 runManager->Initialize();

 
 // Construct the visualization manager
 #ifdef    G4VIS_USE
  G4VisManager * visManager(new G4VisExecutive);
  
  if (visManager==0)
  {
   G4cerr << options.ApplicationName() << ": Vis manager not allocated." << G4endl;
   return;
  }
  
  visManager->Initialize();
 #endif /* G4VIS_USE */
 
 
 // Gets the user interface
 G4UImanager * uiManager(G4UImanager::GetUIpointer());
 
 if (uiManager==0)
 {
  G4cerr << options.ApplicationName() << ": UI manager not allocated." << G4endl;
  return;
 }


 // Construct the messenger to modify the layout
 messenger=new RadmonDetectorMessenger(layout);
 
 if (messenger==0)
 {
  G4cerr << options.ApplicationName() << ": Detector layout messenger not allocated." << G4endl;
  return;
 }
 
 
 // Construct the interactive session
 if (options.Interactive())
 {
  session=new G4UIterminal();
 
  if (session==0)
  {
   G4cerr << options.ApplicationName() << ": Interactive session not allocated." << G4endl;
   return;
  }
 }
 
 
 // Runs startup macros
 RunMacro(options, options.StartupFileName());
 
 if (options.FileName())
  if (!RunMacro(options, options.FileName()))
  {
   G4cerr << options.ApplicationName() << ": File \"" << options.FileName() << "\" not found." << G4endl; 
   return;
  }
  
 // Runs the interactive session
 if (options.Interactive())
 {
  if (options.Verbose())
   G4cout << options.ApplicationName() << ": Interactive session starts ..." << G4endl;
  
  session->SessionStart();
 } 
}



                                                RadmonApplication :: ~RadmonApplication()
{
 // Desctruct the interactive session
 if (session)
  delete session;
  
 // Destruct the messenger to modify the layout 
 if (messenger)
  delete messenger;
 
 // Destruct the visualization manager
 #ifdef    G4VIS_USE
  if (visManager)
   delete visManager;
 #endif /* G4VIS_USE */

 // Destruct the detectors factory
 if (factory)
  delete factory;

 // Destruct the detector layout
 if (layout)
  delete layout;
 
 // Destruct the default run manager
 if (runManager)
  delete runManager;
 
}





G4bool                                          RadmonApplication :: CreateEntityConstructors(void)
{
 return true;
}





G4bool                                          RadmonApplication :: RunMacro(const RadmonApplicationOptions & options, const char * fileName)
{
 std::ifstream test(fileName);
 
 if (!test.good())
  return false;
  
 test.close();
 
 if (options.Verbose())
  G4cout << options.ApplicationName() << ": Running macro \"" << fileName << "\" ..." << G4endl;

 G4String command("/control/execute ");
 command+=fileName;
 uiManager->ApplyCommand(command);
 
 return true;
}
