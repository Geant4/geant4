//
// File name:     RadmonApplication.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplication.cc,v 1.7 2005-10-25 16:39:12 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonApplication.hh"
#include "RadmonApplicationOptions.hh"

#include "RadmonDetectorLayout.hh"
#include "RadmonDetectorConstruction.hh"
#include "RadmonDetectorLabelledEntitiesConstructorsFactory.hh"
#include "RadmonDetectorMessenger.hh"

#include "RadmonPhysicsDummyPhysicsList.hh"

#include "RadmonGeneratorLayout.hh"
#include "RadmonPrimaryGeneratorAction.hh"
#include "RadmonGeneratorsWithLabelFactory.hh"
#include "RadmonGeneratorMessenger.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIdirectory.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef    G4VIS_USE
 #include "G4VisExecutive.hh"
#endif /* G4VIS_USE */



                                                RadmonApplication :: RadmonApplication(const RadmonApplicationOptions & options)
:
 RadmonApplicationDetectorSetup(options),
 RadmonApplicationGeneratorSetup(options),
 valid(false),
 runManager(0),
 detectorLayout(0),
 generatorLayout(0),
 detectorsFactory(0),
 generatorsFactory(0),
 #ifdef    G4VIS_USE
  visManager(0),
 #endif /* G4VIS_USE */
 uiManager(0),
 detectorMessenger(0),
 generatorMessenger(0),
 session(0),
 directory(0)
{
 // Construct the default run manager
 runManager=new G4RunManager();
 
 if (runManager==0)
 {
  G4cerr << options.ApplicationName() << ": Run manager not allocated." << G4endl;
  return;
 }
 
 
 // Construct the detector layout
 detectorLayout=new RadmonDetectorLayout;
 
 if (detectorLayout==0)
 {
  G4cerr << options.ApplicationName() << ": Detector layout not allocated." << G4endl;
  return;
 }
 
 
 // Construct the generator layout
 generatorLayout=new RadmonGeneratorLayout;
 
 if (generatorLayout==0)
 {
  G4cerr << options.ApplicationName() << ": Generator layout not allocated." << G4endl;
  return;
 }
 
 
 // Construct the detectors factory
 detectorsFactory=new RadmonDetectorLabelledEntitiesConstructorsFactory;
 
 if (detectorsFactory==0)
 {
  G4cerr << options.ApplicationName() << ": Detectors factory not allocated." << G4endl;
  return;
 }
 
 
 // Construct the entity constructors
 if (!CreateDetectorEntityConstructors(detectorsFactory))
 {
  G4cerr << options.ApplicationName() << ": Entity constructors not allocated." << G4endl;
  return;
 }
 
 
 // Construct the generators factory
 generatorsFactory=new RadmonGeneratorsWithLabelFactory;
 
 if (generatorsFactory==0)
 {
  G4cerr << options.ApplicationName() << ": Generators factory not allocated." << G4endl;
  return;
 }
 
 
 // Construct the generators algorithms
 if (!CreateGenerators(generatorsFactory))
 {
  G4cerr << options.ApplicationName() << ": Generator algorithms not allocated." << G4endl;
  return;
 }
 
 
 // Construct the detector construction
 RadmonDetectorConstruction * detectorConstruction(new RadmonDetectorConstruction(detectorLayout, detectorsFactory));
 
 if (detectorConstruction==0)
 {
  G4cerr << options.ApplicationName() << ": Detector construction not allocated." << G4endl;
  return;
 }
 
 
 // The detectors factory will is owned by the detectorConstruction
 detectorsFactory=0; 
 
 runManager->SetUserInitialization(detectorConstruction);
 

 // Construct the detector construction
 RadmonPrimaryGeneratorAction * primaryGenerator(new RadmonPrimaryGeneratorAction(generatorLayout, generatorsFactory));
 
 if (detectorConstruction==0)
 {
  G4cerr << options.ApplicationName() << ": Primary generator not allocated." << G4endl;
  return;
 }
 
 
 // The primary generators factory will is owned by the primaryGenerator
 generatorsFactory=0;
 
 runManager->SetUserAction(primaryGenerator);
 
            
 // Set mandatory physics list
 runManager->SetUserInitialization(new RadmonPhysicsDummyPhysicsList);
         

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
 uiManager = G4UImanager::GetUIpointer();
 
 if (uiManager==0)
 {
  G4cerr << options.ApplicationName() << ": UI manager not allocated." << G4endl;
  return;
 }


 // Construct the Radmon directory
 directory = new G4UIdirectory("/radmon/");
 
 if (directory==0)
 {
  G4cerr << options.ApplicationName() << ": Radmon directory not allocated." << G4endl;
  return;
 }
 
 directory->SetGuidance("Radmon application directory.");


 // Construct the messenger to modify the detector layout
 detectorMessenger=new RadmonDetectorMessenger(detectorLayout);
 
 if (detectorMessenger==0)
 {
  G4cerr << options.ApplicationName() << ": Detector layout messenger not allocated." << G4endl;
  return;
 }
 
 
 // Construct the messenger to modify the generator layout
 generatorMessenger=new RadmonGeneratorMessenger(generatorLayout);
 
 if (generatorMessenger==0)
 {
  G4cerr << options.ApplicationName() << ": Generator layout messenger not allocated." << G4endl;
  return;
 }
 
 
 // Construct the interactive session
 if (options.Interactive())
 {
  session=new G4UIterminal(new G4UItcsh);
 
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
 
 valid=true;
}



                                                RadmonApplication :: ~RadmonApplication()
{
 // Desctruct the interactive session
 delete session;
  
 // Destruct the generator messenger to modify the layout 
 delete generatorMessenger;
  
 // Destruct the detector messenger to modify the layout 
 delete detectorMessenger;
  
 // Destruct the Radmon directory
 delete directory;
 
 // Destruct the visualization manager
 #ifdef    G4VIS_USE
  delete visManager;
 #endif /* G4VIS_USE */

 // Destruct the generators factory
 delete generatorsFactory;

 // Destruct the detectors factory
 delete detectorsFactory;

 // Destruct the default run manager
 delete runManager;
 
 // Destruct the generator layout
 delete generatorLayout; 

 // Destruct the detector layout
 delete detectorLayout; 
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
