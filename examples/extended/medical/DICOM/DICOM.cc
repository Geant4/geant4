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
// *******************************************************
#include "G4Types.hh"

#include "G4RunManagerFactory.hh"

#include "globals.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "G4GenericPhysicsList.hh"

#include "DicomRegularDetectorConstruction.hh"
#include "DicomNestedParamDetectorConstruction.hh"
#include "DicomPartialDetectorConstruction.hh"

#include "DicomActionInitialization.hh"

#ifdef G4_DCMTK
#   include "DicomFileMgr.hh"
#else
#   include "DicomHandler.hh"
#endif

#include "DicomIntersectVolume.hh"
#include "QGSP_BIC.hh"
#include "G4tgrMessenger.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Shielding.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Instantiate G4UIExecutive if interactive mode
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  new G4tgrMessenger;
  char* part = std::getenv( "DICOM_PARTIAL_PARAM" );
  G4bool bPartial = FALSE;
  if( part && G4String(part) == "1" ) {
    bPartial = TRUE;
  }

  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(G4long(24534575684783));
  G4long seeds[2];
  seeds[0] = G4long(534524575674523);
  seeds[1] = G4long(526345623452457);
  CLHEP::HepRandom::setTheSeeds(seeds);

  // Construct the default run manager
  char* nthread_c = std::getenv("DICOM_NTHREADS");

  unsigned nthreads = 4;
  unsigned env_threads = 0;

  if(nthread_c) {env_threads=unsigned(G4UIcommand::ConvertToDouble(nthread_c));}
  if(env_threads > 0) {nthreads=env_threads;}

  auto* runManager = G4RunManagerFactory::CreateRunManager();
  runManager->SetNumberOfThreads(nthreads);

  DicomDetectorConstruction* theGeometry = 0;

#ifdef G4_DCMTK
  DicomFileMgr* theFileMgr = 0;
#else
  DicomHandler* dcmHandler = 0;
#endif

  if( !bPartial ){
#ifdef G4_DCMTK

    theFileMgr = DicomFileMgr::GetInstance();
    theFileMgr->Convert("Data.dat");

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
  runManager->SetUserInitialization(phys);

  // Set user action classes
  runManager->SetUserInitialization(new DicomActionInitialization());

  runManager->Initialize();

  new DicomIntersectVolume();

  // visualisation manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (ui)
    {
      UImanager->ApplyCommand("/control/execute vis.mac");
      ui->SessionStart();
      delete ui;
    }
  else
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }

  delete visManager;
  delete runManager;

  if( !bPartial ) {
#ifdef G4_DCMTK
    delete theFileMgr;
#endif
  }

  return 0;
}

