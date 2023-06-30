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
///  \file HadronicNucleusXS.cc
///  \brief Main program,
///         hadronic/FlukaCern/ProcessLevel/CrossSection example.
//
//  Author: G.Hugo, 06 January 2023
//
// -------------------------------------------------------------
//
//      GEANT4 HadronNucleusXS
//
///  This application allows the study of G4 cross-sections, 
///  and in addition, of the FLUKA hadron-nucleus inelastic cross sections.
///
///  The user can printout any particle-material XS. 
///  The XS are exactly the ones defined in any default G4 PhysicsList chosen by the user,
///  or from FLUKA (hadron-nucleus inelastic case).
///
///  In the input file, the user can set:
///  - projectile.
///  - target material (element, compound or even mixture).
///  - plotting options.
///
///  All plots (created via the G4 analysis manager) can be dumped 
///  to any of the usually supported formats (e.g. ROOT format), 
///  but also in a Flair-compatible format.
///
///  NB 1: Unlike the FlukaCern/ProcessLevel/FinalState example,
///  the choice here is to directly use physics lists 
///  (hence under the hood, the processes they define),
///  instead of 'hardcoding' processes of interest.
///  This allows to directly study ALL XS, with no possible discrepancy 
///  with respect to what is defined inside the physics lists.
///
///  NB 2: Note that here, the application is completely independent 
///  from the event loop, gun, detector etc: 
///  the XS printout happend outside of the event loop anyway.
///  Hence, the fakeRun mode is used (setting number of events to 0).
///  This implies that no ActionInitialization is needed 
///  (nor would be used anyway, if ever constructed), 
///  and that the detector is a dummy, empty one.

///  Use: build/HadronNucleusXS all_XS.in FTFP_BERT_HP
///       build/HadronNucleusXS all_XS.in G4_HP_CFLUKAHI
//
// -------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManagerFactory.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"

#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

#ifdef G4_USE_FLUKA
#include "G4_HP_CernFLUKAHadronInelastic_PhysicsList.hh"
#include "FLUKAParticleTable.hh"
#endif

#include "XSHistoManager.hh"

#include "G4Exception.hh"
#include "G4UImanager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int main(G4int argc, char** argv) {

  // Check number of arguments
  if (argc != 3) {
    G4Exception("HadronNucleusXS (main)",
                "Wrong number of input arguments.",
                FatalException,
                "Example use: build/HadronNucleusXS all_XS.in FTFP_BERT_HP");
  }

  // Construct a serial RUN MANAGER.
  std::unique_ptr<G4RunManager> runManager(
     G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly));


  // Empty DETECTOR (compulsory).
  const auto dummyDetector = new DetectorConstruction();
  // The detector is owned by G4RunManager.
  runManager->SetUserInitialization(dummyDetector);


  // Get PHYSICS LIST from command line argument.
  // Default: G4_HP_CFLUKAHI.
  const G4String physicsListName = (argc >= 3 ? argv[2] : "G4_HP_CFLUKAHI");
  const G4bool isFLUKAPhysicsList = (physicsListName == "G4_HP_CFLUKAHI");

  G4VModularPhysicsList* physicsList = nullptr;

  // Create physics list with hadron inelastic processes from FLUKA.
  if (isFLUKAPhysicsList) {
#ifdef G4_USE_FLUKA
    physicsList = new G4_HP_CernFLUKAHadronInelastic_PhysicsList();
#else
    G4Exception("HadronNucleusXS.cc",
      "Wrong compilation mode.",
      FatalException,
      "Requested G4_HP_CernFLUKAHadronInelastic physics list.\n" \
      "This requires COMPILATION in FLUKA mode.\n" \
      "Please fully recompile the example with G4_USE_FLUKA=yes.\n" \
      "For example:\n" \
      "source geant4/examples/extended/hadronic/FlukaCern/FlukaInterface/" \
      "env_FLUKA_G4_interface.sh\n" \
      "cd geant4/examples/extended/hadronic/FlukaCern/ProcessLevel/CrossSection/\n" \
      "mkdir build\n" \
      "cd build\n" \
      "cmake -DGeant4_DIR=your_path_to_geant4 -DG4_USE_FLUKA=1 .. \n" \
      "make -j8 G4_USE_FLUKA=1\n" \
      "NB: First time use of FLUKA interface:\n" \
      "Do not forget to first compile the FLUKA interface itself.\n" \
      "For example: cd geant4/examples/extended/hadronic/FlukaCern/FlukaInterface/ " \
      "&& make interface && make env\n" \
      "FlukaInterface/env_FLUKA_G4_interface.sh then needs to be sourced\n" \
      "in whichever terminal you want to use the FLUKA interface.\n");
#endif
  }
  // Create G4 physics list from available catalog.
  else {
    auto factory = G4PhysListFactory();
    physicsList = factory.GetReferencePhysList(physicsListName);
  }

  // The physics list is owned by G4RunManager.
  runManager->SetUserInitialization(physicsList);

#ifdef G4_USE_FLUKA
  if (isFLUKAPhysicsList) {
    // Initialize FLUKA <-> G4 particles conversions tables.
    fluka_particle_table::initialize();
  }
#endif


  // Create HISTO MANAGER (and its messenger).
  auto histoManager = XSHistoManager();


  // User interface manager (owned by G4RunManagerKernel).
  const auto uiManager = G4UImanager::GetUIpointer();
  // BATCH MODE.
  const G4String command = "/control/execute ";
  const G4String fileName = argv[1];
  uiManager->ApplyCommand(command + fileName);


  // CREATE HISTOGRAMS AND FILL THEM (independent from event loop!).
  histoManager.Book();
  histoManager.EndOfRun();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
