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
// File name:     RadmonApplication.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplication.cc,v 1.12.2.2 2006/06/29 16:08:29 gunter Exp $
// Tag:           $Name: geant4-09-02 $
//

// Include files
#include "RadmonApplication.hh"
#include "RadmonApplicationOptions.hh"

#include "RadmonDetectorLayout.hh"
#include "RadmonDetectorConstruction.hh"
#include "RadmonDetectorLabelledEntitiesConstructorsFactory.hh"
#include "RadmonDetectorMessenger.hh"

#include "RadmonGeneratorLayout.hh"
#include "RadmonPrimaryGeneratorAction.hh"
#include "RadmonGeneratorsWithLabelFactory.hh"
#include "RadmonGeneratorMessenger.hh"

#include "RadmonPhysicsLayout.hh"
#include "RadmonPhysicsList.hh"
#include "RadmonSubPhysicsListWithLabelFactory.hh"
#include "RadmonPhysicsMessenger.hh"

#ifdef G4ANALYSIS_USE
 #include "RadmonAnalysisLayout.hh"
 #include "RadmonAnalysis.hh"
 #include "RadmonDataAnalysisWithLabelFactory.hh"
 #include "RadmonAnalysisMessenger.hh"
#endif /* G4ANALYSIS_USE */

#include "RadmonEventAction.hh"
#include "RadmonSteppingAction.hh"
#include "RadmonApplicationMessenger.hh"

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
 RadmonApplicationPhysicsSetup(options),
 #ifdef    G4ANALYSIS_USE
  RadmonApplicationAnalysisSetup(options),
 #endif /* G4ANALYSIS_USE */
 valid(false),
 runManager(0),
 detectorLayout(0),
 generatorLayout(0),
 physicsLayout(0),
 #ifdef    G4ANALYSIS_USE
  analysisLayout(0),
 #endif /* G4ANALYSIS_USE */
 detectorsFactory(0),
 generatorsFactory(0),
 physicsFactory(0),
 #ifdef    G4ANALYSIS_USE
  analysisFactory(0),
  analysis(0),
 #endif /* G4ANALYSIS_USE */
 #ifdef    G4VIS_USE
  visManager(0),
 #endif /* G4VIS_USE */
 uiManager(0),
 detectorMessenger(0),
 generatorMessenger(0),
 physicsMessenger(0),
 #ifdef    G4ANALYSIS_USE
  analysisMessenger(0),
 #endif /* G4ANALYSIS_USE */
 applicationMessenger(0),
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
 
 
 // Construct the physics list layout
 physicsLayout=new RadmonPhysicsLayout;
 
 if (physicsLayout==0)
 {
  G4cerr << options.ApplicationName() << ": Physics list layout not allocated." << G4endl;
  return;
 }
 
 
 // Construct the analysis list layout
 #ifdef    G4ANALYSIS_USE
  analysisLayout=new RadmonAnalysisLayout;
 
  if (analysisLayout==0)
  {
   G4cerr << options.ApplicationName() << ": Analysis layout not allocated." << G4endl;
   return;
  }
 #endif /* G4ANALYSIS_USE */ 
 
 
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
 
 
 // Construct the physics list factory
 physicsFactory=new RadmonSubPhysicsListWithLabelFactory;
 
 if (physicsFactory==0)
 {
  G4cerr << options.ApplicationName() << ": Physics list factory not allocated." << G4endl;
  return;
 }
 
 
 // Construct the sub physics lists
 if (!CreateSubPhysicsList(physicsFactory))
 {
  G4cerr << options.ApplicationName() << ": Sub physics lists not allocated." << G4endl;
  return;
 }
 
 
 // Construct the analysis list factory
 #ifdef    G4ANALYSIS_USE
  analysisFactory=new RadmonDataAnalysisWithLabelFactory;
 
  if (analysisFactory==0)
  {
   G4cerr << options.ApplicationName() << ": Analysis factory not allocated." << G4endl;
   return;
  }
 
 
  // Construct the data analyses
  if (!CreateDataAnalysis(analysisFactory))
  {
   G4cerr << options.ApplicationName() << ": Data analyses not allocated." << G4endl;
   return;
  }
 #endif /* G4ANALYSIS_USE */
 
 
 // Construct the physics list
 RadmonPhysicsList * physicsList(new RadmonPhysicsList(physicsLayout, physicsFactory));
 
 if (physicsList==0)
 {
  G4cerr << options.ApplicationName() << ": Physics list not allocated." << G4endl;
  return;
 }
 
 
 // The subphysics list factory will be owned by the physicsList
 physicsFactory=0;
 
 runManager->SetUserInitialization(physicsList);
 
            
 // Construct the detector construction
 RadmonDetectorConstruction * detectorConstruction(new RadmonDetectorConstruction(detectorLayout, detectorsFactory));
 
 if (detectorConstruction==0)
 {
  G4cerr << options.ApplicationName() << ": Detector construction not allocated." << G4endl;
  return;
 }
 
 
 // The detectors factory will be owned by the detectorConstruction
 detectorsFactory=0; 
 
 runManager->SetUserInitialization(detectorConstruction);
 

 // Construct the primary generator
 RadmonPrimaryGeneratorAction * primaryGenerator(new RadmonPrimaryGeneratorAction(generatorLayout, generatorsFactory));
 
 if (detectorConstruction==0)
 {
  G4cerr << options.ApplicationName() << ": Primary generator not allocated." << G4endl;
  return;
 }
 
 
 // The primary generators factory will be owned by the primaryGenerator
 generatorsFactory=0;
 
 runManager->SetUserAction(primaryGenerator);
 
            
 // Construct the analysis
 #ifdef    G4ANALYSIS_USE
  AIDA::IAnalysisFactory * aida(AIDA_createAnalysisFactory());
 
  if (aida==0)
  {
   G4cerr << options.ApplicationName() << ": AIDA_createAnalysisFactory returned 0." << G4endl;
   return;
  }
 
  analysis=new RadmonAnalysis(analysisLayout, analysisFactory, aida);
 
  if (analysis==0)
  {
   delete aida;
   G4cerr << options.ApplicationName() << ": Analysis not allocated." << G4endl;
   return;
  }
  
  // The data analyses factory will be owned by the analysis
  analysisFactory=0;
 
  RadmonEventAction * eventAction(RadmonEventAction::Instance());
  eventAction->AttachObserver(analysis);
 #endif /* G4ANALYSIS_USE */
            
 // Initialize the run manager (disabled in order to have UI physics list)
 // runManager->Initialize();

 
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
 
 
 // Construct the messenger to modify the physics list layout
 physicsMessenger=new RadmonPhysicsMessenger(physicsLayout);
 
 if (physicsMessenger==0)
 {
  G4cerr << options.ApplicationName() << ": Physics list layout messenger not allocated." << G4endl;
  return;
 }
 
 
 // Construct the messenger to modify the analysis layout
 #ifdef    G4ANALYSIS_USE
  analysisMessenger=new RadmonAnalysisMessenger(analysisLayout);
 
  if (analysisMessenger==0)
  {
   G4cerr << options.ApplicationName() << ": Analysis layout messenger not allocated." << G4endl;
   return;
  }
 #endif /* G4ANALYSIS_USE */
 
 
 // Construct the messenger to modify applications options
 applicationMessenger=new RadmonApplicationMessenger;
 
 if (applicationMessenger==0)
 {
  G4cerr << options.ApplicationName() << ": Application messenger not allocated." << G4endl;
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

 // Destruct the messenger to modify applications options
 delete applicationMessenger;  

 // Destruct the analysis messenger to modify the layout
 #ifdef    G4ANALYSIS_USE
  delete analysisMessenger;  
 #endif /* G4ANALYSIS_USE */
 
 // Destruct the physics list messenger to modify the layout 
 delete physicsMessenger;
  
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

 // Destruct the analysis
 #ifdef    G4ANALYSIS_USE
  delete analysis;  

  // Destruct the analyses factory
  delete analysisFactory;  
 #endif /* G4ANALYSIS_USE */
 
 // Destruct the sub physics list factory
 delete physicsFactory;

 // Destruct the generators factory
 delete generatorsFactory;

 // Destruct the detectors factory
 delete detectorsFactory;

 // Destruct the default run manager
 delete runManager;
 
 // Destruct the analysis layout
 #ifdef    G4ANALYSIS_USE
  delete analysisLayout; 
 #endif /* G4ANALYSIS_USE */

 // Destruct the physics list layout
 delete physicsLayout; 

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
