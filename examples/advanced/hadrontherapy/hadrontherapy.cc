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
//      ----------------------------------------------------------------------------
//                              GEANT 4 - Hadrontherapy example
//      ----------------------------------------------------------------------------
//
//
//                      ==========>      WEB LINK   <==========
//
//                     http://www.lns.infn.it/link/Hadrontherapy
//
//                      ==========>    MAIN AUTHORS <==========
//
//                       G.A.P. Cirrone(a)*, F.Romano(a), A. Tramontana (a,f)
//
//                      ==========>   PAST AUTHORS  <==========
//
//                      R. Calcagno(a), G.Danielsen (b), F.Di Rosa(a),
//                      S.Guatelli(c), A.Heikkinen(b), P.Kaitaniemi(b),
//                      A.Lechner(d), S.E.Mazzaglia(a),  M.G.Pia(e), G.Russo(a),
//                      M.Russo(a), A.Varisano(a)
//
//
//              (a) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//              (b) Helsinki Institute of Physics, Helsinki, Finland
//              (c) University of Wallongong, Australia
//              (d) CERN, (CH)
//              (e) INFN Section of Genova, genova, Italy
//              (f) Physics and Astronomy Department, Universituy of Catania, Catania, Italy
//
//          *Corresponding author, email to pablo.cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "HadrontherapyEventAction.hh"
#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyDetectorSD.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapyRunAction.hh"
#include "HadrontherapyMatrix.hh"
#include "Randomize.hh"

#include "G4UImessenger.hh"
#include "globals.hh"
#include "HadrontherapySteppingAction.hh"
#include "HadrontherapyAnalysisManager.hh"
#include "HadrontherapyGeometryController.hh"
#include "HadrontherapyGeometryMessenger.hh"
#include "HadrontherapyInteractionParameters.hh"
#include "HadrontherapyLet.hh"

#include "G4ScoringManager.hh"
#include "G4ParallelWorldPhysics.hh"
#include <time.h>

//************************MT*********************
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "HadrontherapyActionInitialization.hh"

//************************MT*********************

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc ,char ** argv)
{
    // Set the Random engine
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
   
  // Only if an initial random seed is needed
    G4int seed =1414159599;// time(0);
  CLHEP::HepRandom::setTheSeed(seed);    
  // G4cout << "******************************************************************"<< seed << G4endl;

    //************************MT*********************
#ifdef G4MULTITHREADED
    
    G4MTRunManager* runManager = new G4MTRunManager;
    //runManager->SetNumberOfThreads(2); // Is equal to 2 by default, it can be setted also with the macro command: /run/numberOfThread 2
#else
    G4RunManager* runManager = new G4RunManager;
#endif
    
    //************************MT*********************
    //   G4RunManager* runManager = new G4RunManager;
    
    // Geometry controller is responsible for instantiating the
    // geometries. All geometry specific setup tasks are now in class
    // HadrontherapyGeometryController.
    HadrontherapyGeometryController *geometryController = new HadrontherapyGeometryController();
	
    // Connect the geometry controller to the G4 user interface
    HadrontherapyGeometryMessenger *geometryMessenger = new HadrontherapyGeometryMessenger(geometryController);
	
    G4ScoringManager *scoringManager = G4ScoringManager::GetScoringManager();
    scoringManager->SetVerboseLevel(1);
    
	
    // Initialize the default Hadrontherapy geometry
    geometryController->SetGeometry("default");
	
    // Initialize command based scoring
    G4ScoringManager::GetScoringManager();
	
    // Initialize the physics
    G4PhysListFactory factory;
    G4VModularPhysicsList* phys = 0;
    G4String physName = "";
    
    // Physics List name defined via environment variable
    char* path = getenv("PHYSLIST");
    if (path) { physName = G4String(path); }
    
    if(physName != "" && factory.IsReferencePhysList(physName))
    {
        phys = factory.GetReferencePhysList(physName);
    }
    if (phys)
    {
        G4cout << "Going to register G4ParallelWorldPhysics" << G4endl;
        phys->RegisterPhysics(new G4ParallelWorldPhysics("DetectorROGeometry"));
    }
    else
    {
        G4cout << "Using HadrontherapyPhysicsList()" << G4endl;
        phys = new HadrontherapyPhysicsList();
    }
    
    runManager->SetUserInitialization(phys);
    
    //************************MT
    runManager->SetUserInitialization(new HadrontherapyActionInitialization);
    //************************MT
    
    
    //************************MT: DA SPOSTARE IN hADRONTHERAPYACTIONiNITIALIZATION.CC*********************
    /*
     // Initialize the primary particles
     HadrontherapyPrimaryGeneratorAction *pPrimaryGenerator = new HadrontherapyPrimaryGeneratorAction();
     runManager -> SetUserAction(pPrimaryGenerator);
     
     // Optional UserActions: run, event, stepping
     HadrontherapyRunAction* pRunAction = new HadrontherapyRunAction();
     runManager -> SetUserAction(pRunAction);
     
     HadrontherapyEventAction* pEventAction = new HadrontherapyEventAction();
     runManager -> SetUserAction(pEventAction);
     
     HadrontherapySteppingAction* steppingAction = new HadrontherapySteppingAction(pRunAction);
     runManager -> SetUserAction(steppingAction);
     */
    // Interaction data: stopping powers
    HadrontherapyInteractionParameters* pInteraction = new HadrontherapyInteractionParameters(true);
	
    // Initialize analysis
    HadrontherapyAnalysisManager* analysis = HadrontherapyAnalysisManager::GetInstance();
#ifdef G4ANALYSIS_USE_ROOT
    analysis -> book();
#endif
    
    // Get the pointer to the visualization manager
#ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
    visManager -> Initialize();
#endif
    
    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    
    // If no macro file is passed as argument,
    // the User Interface is called
    if (argc==1)
    {
#ifdef G4UI_USE
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
        
        if(factory.IsReferencePhysList(physName))
        {
            UImanager->ApplyCommand("/control/execute macro/defaultMacroWithReferencePhysicsList.mac");
        }
        else
        {
            UImanager->ApplyCommand("/control/execute macro/defaultMacro.mac");
        }
        
#endif
        ui->SessionStart();
        delete ui;
#endif
    }
    // Batch mode: the following commands are called when a macro file
    // is passed as argument
    else
    {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
    }
    
    // Job termination
    // Store dose & fluence data to ASCII & ROOT files
    if ( HadrontherapyMatrix * pMatrix = HadrontherapyMatrix::GetInstance() )
    {
        pMatrix -> TotalEnergyDeposit();
        pMatrix -> StoreDoseFluenceAscii();
#ifdef G4ANALYSIS_USE_ROOT
        
        pMatrix -> StoreDoseFluenceRoot();
#endif
    }
    
    if (HadrontherapyLet *let = HadrontherapyLet::GetInstance())
        if(let -> doCalculation)
        {
            let -> LetOutput(); 	// Calculate let
            let -> StoreLetAscii(); // Store it
#ifdef G4ANALYSIS_USE_ROOT
            
            let -> StoreLetRoot();
#endif
        }
    
    
#ifdef G4ANALYSIS_USE_ROOT
    if (analysis -> IsTheTFile()) analysis -> flush();     // Finalize & write the root file
#endif
    
    
#ifdef G4VIS_USE
    delete visManager;
#endif
    
    delete geometryMessenger;
    delete geometryController;
    delete pInteraction;
    delete runManager;
    delete analysis;
    return 0;
    
}
