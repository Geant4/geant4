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
<<<<<<< HEAD
=======
//                                      MAIN AUTHOR
//                                  ====================
//                                  G.A.P. Cirrone(a)*
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
//                      ==========>      WEB LINK   <==========
//
<<<<<<< HEAD
//                     http://www.lns.infn.it/link/Hadrontherapy
=======
//                                  ACTUAL CONTRIBUTORS
//                                  ====================
//                  G.A.P. Cirrone(a), L. Pandola(a), G. Petringa(a)
//
//                      ==========>    MAIN AUTHORS <==========
//
//                       G.A.P. Cirrone(a)*, F.Romano(a), A. Tramontana (a,f)
//
//                      ==========>   PAST AUTHORS  <==========
//
//                      R. Calcagno(a), G.Danielsen (b), F.Di Rosa(a),
//                      S.Guatelli(c), A.Heikkinen(b), P.Kaitaniemi(b),
//                      A.Lechner(d), S.E.Mazzaglia(a), Z. Mei(h), M.G.Pia(e),
//                      F.Romano(a), G.Russo(a,g), M.Russo(a), A. Tramontana (a),
//                      A.Varisano(a)
//
//              (a) Laboratori Nazionali del Sud of INFN, Catania, Italy
//              (b) Helsinki Institute of Physics, Helsinki, Finland
//              (c) University of Wallongong, Australia
//              (d) CERN, Geneve, Switzwerland
//              (e) INFN Section of Genova, Genova, Italy
//              (f) Physics and Astronomy Department, Univ. of Catania, Catania, Italy
//              (g) CNR-IBFM, Italy
//              (h) Institute of Applied Electromagnetic Engineering(IAEE)
//                  Huazhong University of Science and Technology(HUST), Wuhan, China
//
//
//                                          WEB
//                                      ===========
//       https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy
//
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
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
#include "G4Timer.hh"
#include "G4RunManagerFactory.hh"
#include "HadrontherapyActionInitialization.hh"

<<<<<<< HEAD
//************************MT*********************

#ifdef G4VIS_USE
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc ,char ** argv)
{
<<<<<<< HEAD
=======
        G4UIExecutive* ui = 0;
    if ( argc == 1 ) {
        ui = new G4UIExecutive(argc, argv);
    }
    
    //Instantiate the G4Timer object, to monitor the CPU time spent for
    //the entire execution
    G4Timer* theTimer = new G4Timer();
    //Start the benchmark
    theTimer->Start();
    
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
    // Set the Random engine
    // The following guarantees random generation also for different runs
    // in multithread
    CLHEP::RanluxEngine defaultEngine( 1234567, 4 );
    G4Random::setTheEngine( &defaultEngine );
    G4int seed = (G4int) time( NULL );
    G4Random::setTheSeed( seed );
 
 auto* runManager = G4RunManagerFactory::CreateRunManager();
 G4int nThreads = 4;
 runManager->SetNumberOfThreads(nThreads); 

    // Geometry controller is responsible for instantiating the
    // geometries. All geometry specific m tasks are now in class
    // HadrontherapyGeometryController.
    HadrontherapyGeometryController *geometryController = new HadrontherapyGeometryController();
	
    // Connect the geometry controller to the G4 user interface
    HadrontherapyGeometryMessenger *geometryMessenger = new HadrontherapyGeometryMessenger(geometryController);
	
    G4ScoringManager *scoringManager = G4ScoringManager::GetScoringManager();
    scoringManager->SetVerboseLevel(1);
    
<<<<<<< HEAD
	
    // Initialize the default Hadrontherapy geometry
    geometryController->SetGeometry("default");
	
    // Initialize command based scoring
    G4ScoringManager::GetScoringManager();
	
=======
    // Initialize the default Hadrontherapy geometry
    geometryController->SetGeometry("default");
    
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
    // Initialize the physics
    G4PhysListFactory factory;
    G4VModularPhysicsList* phys = 0;
    G4String physName = "";
    
    // Physics List name defined via environment variable
    char* path = std::getenv("PHYSLIST");
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
    
    // Initialisations of physics
    runManager->SetUserInitialization(phys);
    
    // Initialisation of the Actions
    runManager->SetUserInitialization(new HadrontherapyActionInitialization);
    //************************MT
    
<<<<<<< HEAD
    
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
=======
    // Initialize command based scoring
    G4ScoringManager::GetScoringManager();
    
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
    // Interaction data: stopping powers
    HadrontherapyInteractionParameters* pInteraction = new HadrontherapyInteractionParameters(true);
	
    // Initialize analysis
<<<<<<< HEAD
    HadrontherapyAnalysisManager* analysis = HadrontherapyAnalysisManager::GetInstance();
#ifdef G4ANALYSIS_USE_ROOT
    analysis -> book();
#endif
=======
    HadrontherapyAnalysis* analysis = HadrontherapyAnalysis::GetInstance();
    

>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
    
// Initialise the Visualisation
    G4VisManager* visManager = new G4VisExecutive;
    visManager -> Initialize();
    
    //** Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    
<<<<<<< HEAD
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
=======
    if ( !ui ) {
        // batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
        
    }
    
    else {

        UImanager -> ApplyCommand("/control/execute macro/defaultMacro.mac");
        ui -> SessionStart();
        delete ui;
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
    }
    delete visManager;
 
    //Stop the benchmark here
    theTimer->Stop();
    
    G4cout << "The simulation took: " << theTimer->GetRealElapsed() << " s to run (real time)"
    << G4endl;
    
    // Job termination
<<<<<<< HEAD
    // Store dose & fluence data to ASCII & ROOT files
    if ( HadrontherapyMatrix * pMatrix = HadrontherapyMatrix::GetInstance() )
=======
    // Free the store: user actions, physics_list and detector_description are
    // owned and deleted by the run manager, so they should not be deleted
    // in the main() program !
    
    
        if ( HadrontherapyMatrix * pMatrix = HadrontherapyMatrix::GetInstance() )
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
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
    
<<<<<<< HEAD
    
#ifdef G4ANALYSIS_USE_ROOT
    if (analysis -> IsTheTFile()) analysis -> flush();     // Finalize & write the root file
#endif
    
    
#ifdef G4VIS_USE
    delete visManager;
#endif
    
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
    delete geometryMessenger;
    delete geometryController;
    delete pInteraction;
    delete runManager;
    delete analysis;
    return 0;
    
}
