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
/// \file dicom2.cc
/// \brief Main program of the Dicom2 example

#include "G4Types.hh"
#ifdef G4MULTITHREADED
#   include "G4MTRunManager.hh"
#else
#   include "G4RunManager.hh"
#endif
#include "G4UImanager.hh"
#include "G4GenericPhysicsList.hh"
#include "G4tgrMessenger.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4Timer.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "Shielding.hh"
#include "QGSP_BIC.hh"
#include "QBBC.hh"

#include "DicomRegularDetectorConstruction.hh"
#include "DicomNestedParamDetectorConstruction.hh"
#include "DicomPartialDetectorConstruction.hh"
#include "Dicom2ActionInitialization.hh"
#include "DicomIntersectVolume.hh"
#ifdef G4_DCMTK
#   include "DicomFileMgr.hh"
#else
#   include "DicomHandler.hh"
#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
    // Detect interactive mode (if no arguments) and define UI session
    //
    G4UIExecutive* ui = 0;
    if (argc == 1)
        ui = new G4UIExecutive(argc, argv);

    new G4tgrMessenger;
    char* part = std::getenv("DICOM_PARTIAL_PARAM");
    G4bool bPartial = (part && G4String(part) == "1") ? true : false;

    CLHEP::HepRandom::setTheEngine(new CLHEP::MixMaxRng);
    CLHEP::HepRandom::setTheSeed(G4long(24534575684783));
    G4long seeds[2];
    seeds[0] = G4long(534524575674523);
    seeds[1] = G4long(526345623452457);
    CLHEP::HepRandom::setTheSeeds(seeds);

    // Construct the default run manager
#ifdef G4MULTITHREADED
    G4int nthreads = G4GetEnv<G4int>("DICOM_NTHREADS", G4Thread::hardware_concurrency());
    G4MTRunManager* runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(nthreads);

    G4cout << "\n\n\tDICOM2 running in multithreaded mode with "
           << runManager->GetNumberOfThreads()
           << " threads\n\n" << G4endl;
#else
    G4RunManager* runManager = new G4RunManager;
    G4cout << "\n\n\tDICOM running in serial mode\n\n" << G4endl;
#endif

    DicomDetectorConstruction* theGeometry = 0;

#ifdef G4_DCMTK
    DicomFileMgr* theFileMgr = 0;
#else
    DicomHandler* dcmHandler = 0;
#endif

    if( !bPartial )
    {
#ifdef G4_DCMTK
        G4String inpfile = "Data.dat";
        char* env_inpfile = std::getenv("DICOM_INPUT_FILE");
        if(env_inpfile)
            inpfile = env_inpfile;

        theFileMgr = DicomFileMgr::GetInstance();
        theFileMgr->Convert(inpfile);
#else
        // Treatment of DICOM images before creating the G4runManager
        dcmHandler = DicomHandler::Instance();
        dcmHandler->CheckFileFormat();
#endif

        // Initialisation of physics, geometry, primary particles ...
        char* nest = std::getenv( "DICOM_NESTED_PARAM" );
        if( nest && G4String(nest) == "1" ) {
            theGeometry = new DicomNestedParamDetectorConstruction();
        } else {
            theGeometry = new DicomRegularDetectorConstruction();
        }
    } else {
        theGeometry = new DicomPartialDetectorConstruction();
    }
    runManager->SetUserInitialization(theGeometry);

    //    std::vector<G4String>* MyConstr = new std::vector<G4String>;
    //    MyConstr->push_back("G4EmStandardPhysics");
    //    G4VModularPhysicsList* phys = new G4GenericPhysicsList(MyConstr);
    G4VModularPhysicsList* phys = new Shielding();
    phys->SetDefaultCutValue(0.5*CLHEP::mm);
    runManager->SetUserInitialization(phys);

    // User action initialization
    runManager->SetUserInitialization(new Dicom2ActionInitialization());

    runManager->Initialize();

    new DicomIntersectVolume();

    // Initialize visualization
    //
    G4VisManager* visManager = new G4VisExecutive;
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    G4Timer t;
    t.Start();

    // Process macro or start UI session
    //
    if ( ! ui ) {
        // batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
    }
    else {
        // interactive mode
        UImanager->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
        delete ui;
    }

    t.Stop();

    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    // owned and deleted by the run manager, so they should not be deleted
    // in the main() program !

    delete visManager;
    delete runManager;

    if( !bPartial )
    {
#ifdef G4_DCMTK
        delete theFileMgr;
#endif
    }

    G4cout << "\n[" << argv[0] << "] Primary execution time: " << t << "\n"
           << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
