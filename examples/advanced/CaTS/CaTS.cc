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
//---------------------------------------------------------------------
//*            |\___/|                                                *
//*            )     (                                                *
//*           =\     /=                                               *
//*             )===(                                                 *
//*            /     \     CaTS: Calorimeter and Tracker Simulation   *
//*            |     |     is a flexible and extend-able framework    *
//*           /       \    for the simulation of various detector     *
//*	      \       /    systems                                    *
//*            \__  _/     https://github.com/hanswenzel/CaTS         *
//*	         ( (                                                  *
//*	          ) )                                                 *
//*              (_(                                                  *
//* CaTS also serves as an example that demonstrates how to use       *
//* opticks from within Geant4 for the creation and propagation of    *
//* optical photons.                                                  *
//* see https://bitbucket.org/simoncblyth/opticks.git).               *
//* Ascii Art by Joan Stark: https://www.asciiworld.com/-Cats-2-.html *
//---------------------------------------------------------------------
//
// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file CaTS.cc
/// \brief main driver of CaTS
//
// project headers:
#include "ActionInitialization.hh"
#include "CaTSVersion.hh"
#include "ConfigurationManager.hh"
#include "DetectorConstruction.hh"
#include "PhysicsConfigurator.hh"
// Geant4 headers:
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4Timer.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VisExecutive.hh"
#include <G4Threading.hh>
#ifdef WITH_G4OPTICKS
#  include "OPTICKS_LOG.hh"
#endif
#include "TROOT.h"
#include <thread>

int main(int argc, char** argv)
{
#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif
  G4bool interactive   = false;
  G4String physicsconf = "";
  G4String gdmlfile    = "";
  G4String macrofile   = "";
  G4UIExecutive* ui    = nullptr;
  for(G4int i = 1; i < argc; i = i + 2)
  {
    if(G4String(argv[i]) == "-g")
    {
      gdmlfile = argv[i + 1];
    }
    else if(G4String(argv[i]) == "-pl")
    {
      physicsconf = G4String(argv[i + 1]);
    }
    else if(G4String(argv[i]) == "-m")
    {
      macrofile = G4String(argv[i + 1]);
    }
#ifdef G4MULTITHREADED
    else if(G4String(argv[i]) == "-t")
    {
      nThreads = G4UIcommand::ConvertToInt(argv[i + 1]);
    }
#endif
  }
  if(gdmlfile == "")
  {
    G4cout << "Error! Mandatory input file is not specified!" << G4endl;
    G4cout << G4endl;
    G4cout << G4endl;
    G4cout << "Usage:  CaTS -g input_gdml_file:mandatory" << G4endl;
    G4cout << G4endl;
    return -1;
  }
  G4cout
    << G4endl
    << "---------------------------------------------------------------------"
    << G4endl
    << "*            |\\___/|                                                "
       "*"
    << G4endl
    << "*            )     (                                                *"
    << G4endl
    << "*           =\\     /=                                               "
       "*"
    << G4endl
    << "*             )===(      Welcome to:                                *"
    << G4endl
    << "*            /     \\     CaTS: Calorimeter and Tracker Simulation   "
       "*"
    << G4endl
    << "*            |     |     a flexible and extend-able framework       *"
    << G4endl
    << "*           /       \\    for the simulation of various detector     "
       "*"
    << G4endl
    << "*	    \\       /    systems                                    *"
    << G4endl
    << "*            \\__  _/     https://github.com/hanswenzel/CaTS         "
       "*"
    << G4endl
    << "*              ( (                                                  *"
    << G4endl << "*	        ) )      Version: " << CaTSVersion
    << "                          *" << G4endl
    << "*              (_(       Date:    " << CaTSDate << "                 *"
    << G4endl
    << "---------------------------------------------------------------------"
    << G4endl << G4endl;
  if(physicsconf == "")
  {
    G4cout << "Warning! no physics configuration specified!" << G4endl;
    G4cout << "Using default FTFP_BERT+OPTICAL+STEPLIMIT" << G4endl;
    physicsconf = "FTFP_BERT+OPTICAL+STEPLIMIT";
    G4cout << "Usage:  CaTS -pl physicsconfiguration" << G4endl;
    G4cout << G4endl;
  }
  if(macrofile == "")
  {
    G4cout << "Warning! no macro specified!" << G4endl;
    G4cout << "assume interactive mode" << G4endl;
    interactive = true;
    ui          = new G4UIExecutive(argc, argv);
    G4cout << G4endl;
    G4cout << G4endl;
    G4cout << "Usage:  CaTS -m macrofile" << G4endl;
    G4cout << G4endl;
  }
  G4Timer* eventTimer = new G4Timer;
  eventTimer->Start();
#ifdef WITH_G4OPTICKS
  OPTICKS_LOG(argc, argv);
#endif
  G4VModularPhysicsList* phys =
    PhysicsConfigurator::getInstance()->Construct(physicsconf);
  G4String DumpFilename = gdmlfile + "_G4";
  ConfigurationManager::getInstance()->setGDMLFileName(DumpFilename);
  DetectorConstruction* dc = new DetectorConstruction(gdmlfile);
  // Run manager
  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);
#ifdef G4MULTITHREADED
  // number of threads not set so use number of cores
  if(nThreads == 0)
  {
    nThreads = G4Threading::G4GetNumberOfCores();
  }
  if(nThreads > 0)
  {
    rm->SetNumberOfThreads(nThreads);
  }
#endif
  rm->SetUserInitialization(dc);
  rm->SetUserInitialization(phys);
  ActionInitialization* actionInitialization = new ActionInitialization();
  rm->SetUserInitialization(actionInitialization);
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  if(interactive)
  {
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
    delete visManager;
  }
  else
  {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + macrofile);
    
    delete ui;
  }
  eventTimer->Stop();
  double totalCPUTime =
    eventTimer->GetUserElapsed() + eventTimer->GetSystemElapsed();
  G4int precision_t          = G4cout.precision(3);
  std::ios::fmtflags flags_t = G4cout.flags();
  G4cout.setf(std::ios::fixed, std::ios::floatfield);
  G4cout << "TimeTotal> " << eventTimer->GetRealElapsed() << " " << totalCPUTime
         << G4endl;
  G4cout.setf(flags_t);
  G4cout.precision(precision_t);
  delete eventTimer;
  delete rm;
  return 0;
}
