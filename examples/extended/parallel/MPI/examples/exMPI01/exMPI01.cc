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
/// @file exMPI01.cc
//
/// @brief A MPI example code

#include "G4MPImanager.hh"
#include "G4MPIsession.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"

#include "G4VisExecutive.hh"

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "FTFP_BERT.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // random engine
  CLHEP::MTwistEngine randomEngine;
  G4Random::setTheEngine(&randomEngine);

  // --------------------------------------------------------------------
  // MPI session
  // --------------------------------------------------------------------
  // At first, G4MPImanager/G4MPIsession should be created.
  G4MPImanager* g4MPI = new G4MPImanager(argc, argv);

  // MPI session (G4MPIsession) instead of G4UIterminal
  // Terminal availability depends on your MPI implementation.
  G4MPIsession* session = g4MPI-> GetMPIsession();

  // LAM/MPI users can use G4tcsh.
  G4String prompt = "[40;01;33m";
  prompt += "G4MPI";
  prompt += "[40;31m(%s)[40;36m[%/][00;30m:";
  session-> SetPrompt(prompt);

  // --------------------------------------------------------------------
  // user application setting
  // --------------------------------------------------------------------
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager();
  runManager-> SetNumberOfThreads(4);
#else
  G4RunManager* runManager = new G4RunManager();
#endif

  // setup your application
  runManager-> SetUserInitialization(new DetectorConstruction);
  runManager-> SetUserInitialization(new FTFP_BERT);
  runManager-> SetUserInitialization(new ActionInitialization);

  runManager-> Initialize();

  G4VisExecutive* visManager = new G4VisExecutive;
  visManager-> Initialize();
  G4cout << G4endl;

  // --------------------------------------------------------------------
  // ready for go
  // MPIsession treats both interactive and batch modes.
  // Just start your session as below.
  // --------------------------------------------------------------------
  session-> SessionStart();

  // --------------------------------------------------------------------
  // termination
  // --------------------------------------------------------------------
  delete visManager;

  delete g4MPI;

  delete runManager;

  return EXIT_SUCCESS;
}
