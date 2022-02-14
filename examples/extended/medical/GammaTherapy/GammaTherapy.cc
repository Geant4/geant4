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
<<<<<<< HEAD
// $Id: GammaTherapy.cc 67994 2013-03-13 11:05:39Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
/// \file medical/GammaTherapy/GammaTherapy.cc
/// \brief Main program of the medical/GammaTherapy example
//
//
// -------------------------------------------------------------
//      GEANT4 ibrem
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- ibrem ----------------------------
//
//  Modified: 05.04.01 Vladimir Ivanchenko new design of ibrem
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4Types.hh"

#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "TrackingAction.hh"
#include "RunAction.hh"
#include "Histo.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int main(int argc,char** argv) {

  //detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if (argc == 1) ui = new G4UIExecutive(argc,argv);

  // Construct the default run manager
  auto* runManager = G4RunManagerFactory::CreateRunManager();
  G4int nThreads = std::min(G4Threading::G4GetNumberOfCores(),2);
  if (argc==3) nThreads = G4UIcommand::ConvertToInt(argv[2]);
  runManager->SetNumberOfThreads(nThreads);
  G4cout << "===== GammaTherapy is started with "
         <<  runManager->GetNumberOfThreads() << " threads =====" << G4endl;

  // set mandatory initialization classes
  DetectorConstruction* det = new DetectorConstruction();
  runManager->SetUserInitialization(det);

  runManager->SetUserInitialization(new PhysicsList());

  // set user action classes
  runManager->SetUserAction(new PrimaryGeneratorAction(det));
  runManager->SetUserAction(new RunAction());
  runManager->SetUserAction(new TrackingAction());

  // get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  if (argc==1)   // Define UI terminal for interactive mode
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
  else if (argc>1) // Batch mode with 1 or more files
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }

#ifdef G4VIS_USE
  delete visManager;
#endif

  // job termination
  delete runManager;

  return 0;
}


