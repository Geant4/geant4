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
// $Id: DICOM.cc 74809 2013-10-22 09:49:26Z gcosmo $
//
/// \file medical/DICOM/DICOM.cc
/// \brief Main program of the medical/DICOM example
//
// The code was written by :
//        *Louis Archambault louis.archambault@phy.ulaval.ca,
//      *Luc Beaulieu beaulieu@phy.ulaval.ca
//      +Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
//
// *Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//
// + Université Laval, Québec (QC) Canada
//*******************************************************//


#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "globals.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DicomPhysicsList.hh"

#include "DicomRegularDetectorConstruction.hh"
#include "DicomNestedParamDetectorConstruction.hh"
#include "DicomPartialDetectorConstruction.hh"

#include "DicomActionInitialization.hh"

#include "DicomHandler.hh"
#include "DicomIntersectVolume.hh"
#include "QGSP_BIC.hh"
#include "G4tgrMessenger.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//=================================================================================

int main(int argc,char** argv)
{

    new G4tgrMessenger;
    char* part = getenv( "DICOM_PARTIAL_PARAM" );
    G4bool bPartial = FALSE;
    if( part && G4String(part) == "1" ) {
        bPartial = TRUE;
    }

    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheSeed(24534575684783);
    long seeds[2];
    seeds[0] = 534524575674523;
    seeds[1] = 526345623452457;
    CLHEP::HepRandom::setTheSeeds(seeds);

    // Construct the default run manager
#ifdef G4MULTITHREADED
    char* nthread_c = getenv("DICOM_NTHREADS");

    unsigned nthreads = 4;
    unsigned env_threads = 0;
    
    if(nthread_c) { env_threads = G4UIcommand::ConvertToDouble(nthread_c); }
    if(env_threads > 0) { nthreads = env_threads; }

    G4MTRunManager* runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(nthreads);

    std::cout << "\n\n\tDICOM running in multithreaded mode with " << nthreads 
          << " threads\n\n" << std::endl;

    
#else
    G4RunManager* runManager = new G4RunManager;
    std::cout << "\n\n\tDICOM running in serial mode\n\n" << std::endl;

#endif
    
    DicomDetectorConstruction* theGeometry = 0;
    DicomHandler* dcmHandler = 0;

    if( !bPartial ){
        // Treatment of DICOM images before creating the G4runManager
        dcmHandler = new DicomHandler;
        dcmHandler->CheckFileFormat();

        // Initialisation of physics, geometry, primary particles ...
        char* nest = getenv( "DICOM_NESTED_PARAM" );
        if( nest && G4String(nest) == "1" ) {
            theGeometry = new DicomNestedParamDetectorConstruction();
        } else {
            theGeometry = new DicomRegularDetectorConstruction();
        }
    } else {
        theGeometry = new DicomPartialDetectorConstruction();
    }

    // runManager->SetUserInitialization(new DicomPhysicsList);
    G4VUserPhysicsList* phys = 0;
    runManager->SetUserInitialization(phys = new QGSP_BIC);
    phys->SetVerboseLevel(0);
    runManager->SetUserInitialization(theGeometry);

    // Set user action classes
    runManager->SetUserInitialization(new DicomActionInitialization());

    runManager->Initialize();

    new DicomIntersectVolume();

#ifdef G4VIS_USE
    // visualisation manager
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
#endif


    G4UImanager* UImanager = G4UImanager::GetUIpointer();


    if (argc==1)
    {
#ifdef G4UI_USE
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
        UImanager->ApplyCommand("/control/execute vis.mac");
#endif
        ui->SessionStart();
        delete ui;
#endif
    }
    else
    {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
    }
    
    delete runManager;
    
#ifdef G4VIS_USE
    delete visManager;
#endif
    
    if( !bPartial ) { delete dcmHandler; }
    
    return 0;
}



