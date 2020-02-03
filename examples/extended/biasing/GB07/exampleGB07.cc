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
/// \file GB07/exampleGB07.cc
/// \brief Main program of the GB07 example
//
//

#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "GB07ActionInitialization.hh"

#include "G4UImanager.hh"

#include "GB07DetectorConstruction.hh"
#include "GB07PrimaryGeneratorAction.hh"

#include "FTFP_BERT.hh"
#include "G4GenericBiasingPhysics.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " ./exampleGB07 [-m macro ] "
           << " [-b biasing {'on','off'}]"
           << "\n or\n ./exampleGB07 [macro.mac]"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  if ( argc > 5 ) {
    PrintUsage();
    return 1;
  }

  G4String macro("");
  G4String onOffBiasing("");
  if ( argc == 2 ) macro = argv[1];
  else
    {
      for ( G4int i=1; i<argc; i=i+2 )
        {
          if      ( G4String(argv[i]) == "-m" ) macro        = argv[i+1];
          else if ( G4String(argv[i]) == "-b" ) onOffBiasing = argv[i+1];
          else
            {
              PrintUsage();
              return 1;
            }
        }
    }

  if ( onOffBiasing == "" ) onOffBiasing = "on";

  // Instantiate G4UIExecutive if interactive mode
  G4UIExecutive* ui = nullptr;
  if ( macro == "" ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // -- Construct the run manager : MT or sequential one
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  G4cout << "      ********** Run Manager constructed in MT mode ************ " << G4endl;
  // -- Choose 4 threads:
  runManager->SetNumberOfThreads(4);
#else
  G4RunManager * runManager = new G4RunManager;
  G4cout << "      ********** Run Manager constructed in sequential mode ************ " << G4endl;
#endif

  // -- Set mandatory initialization classes
  GB07DetectorConstruction* detector = new GB07DetectorConstruction();
  runManager->SetUserInitialization(detector);
  // -- Select a physics list
  FTFP_BERT* physicsList = new FTFP_BERT;
  // -- And augment it with biasing facilities:
  G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics();
  if ( onOffBiasing == "on" )
    {
      // -- Specify the processes that will be under biasing:
      // ---- Decide to apply the technique to hadronic inelastic interactions:
      std::vector< G4String >
         piPlusProcessesToBias,     piMinusProcessesToBias,
   	 protonProcessesToBias,  antiProtonProcessesToBias,
	neutronProcessesToBias, antiNeutronProcessesToBias;
      piPlusProcessesToBias     .push_back(         "pi+Inelastic");
      piMinusProcessesToBias    .push_back(         "pi-Inelastic");
      protonProcessesToBias     .push_back(      "protonInelastic");
      antiProtonProcessesToBias .push_back( "anti_protonInelastic");
      neutronProcessesToBias    .push_back(     "neutronInelastic");
      neutronProcessesToBias    .push_back("nCapture"             );
      antiNeutronProcessesToBias.push_back("anti_neutronInelastic");
      biasingPhysics->PhysicsBias("pi+",               piPlusProcessesToBias);
      biasingPhysics->PhysicsBias("pi-",              piMinusProcessesToBias);
      biasingPhysics->PhysicsBias("proton",            protonProcessesToBias);
      biasingPhysics->PhysicsBias("anti_proton",   antiProtonProcessesToBias);
      biasingPhysics->PhysicsBias("neutron",          neutronProcessesToBias);
      biasingPhysics->PhysicsBias("anti_neutron", antiNeutronProcessesToBias);
      // ---- Apply also the technique to some EM processes (producing at least
      // ---- 2 same or 2 particle - anti-particle pairs):
      std::vector< G4String > gammaProcessesToBias,
	electronProcessesToBias, positronProcessesToBias;
      gammaProcessesToBias      .push_back( "conv"            );
      gammaProcessesToBias      .push_back( "photonNuclear"   );
      electronProcessesToBias   .push_back( "electronNuclear" );
      positronProcessesToBias   .push_back( "annihil"         );
      positronProcessesToBias   .push_back( "positronNuclear" );
      biasingPhysics->PhysicsBias("gamma",              gammaProcessesToBias);
      biasingPhysics->PhysicsBias("e-",              electronProcessesToBias);
      biasingPhysics->PhysicsBias("e+",              positronProcessesToBias);
      // ---- And apply it to pi0 decay:
      std::vector< G4String > pi0ProcessesToBias;
      pi0ProcessesToBias        .push_back( "Decay"           );
      biasingPhysics->PhysicsBias("pi0",                  pi0ProcessesToBias);
      // --
      physicsList->RegisterPhysics(biasingPhysics);
      G4cout << "      ********************************************************* " << G4endl;
      G4cout << "      ********** processes are wrapped for biasing ************ " << G4endl;
      G4cout << "      ********************************************************* " << G4endl;
    }
  else
    {
      G4cout << "      ************************************************* " << G4endl;
      G4cout << "      ********** processes are not wrapped ************ " << G4endl;
      G4cout << "      ************************************************* " << G4endl;
    }
  runManager->SetUserInitialization(physicsList);
  // -- Action initialization:
  runManager->SetUserInitialization(new GB07ActionInitialization);

  // Initialize G4 kernel
  runManager->Initialize();

  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( macro != "" )   // batch mode
    {
      G4String command = "/control/execute ";
      UImanager->ApplyCommand(command+macro);
    }
  else
    {  // interactive mode : define UI session
      UImanager->ApplyCommand("/control/execute vis.mac");
      //      if (ui->IsGUI())
      //        UImanager->ApplyCommand("/control/execute gui.mac");
      ui->SessionStart();
      delete ui;
    }

  delete visManager;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
