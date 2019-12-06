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
/// \file extensibleFactory.cc
/// \brief Main program of the extensibleFactory example
//
//
//
// -------------------------------------------------------------
//      extensibleFactory
//
//  Application demonstrating the extensible physics list factory
//
//  Author of hadronic/Hadr00/Hadr00.cc
//     V.Ivanchenko, 20 June 2008  (as hadronic/Hadr00/Hadr00.cc)
//  Author of examples/extended/physicslists/factory/factory.cc
//      I. Hrivnacova, 2017-09-26

//  Modified from factory.cc
//      R.Hatcher 2017-10-31
//        copied from examples/extended/physicslists/factory
//        modified to use alternative extensible physics list factory
//
// -------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

/////////////////////////////////////////////////////////////////////////////
// The following change is the _only_ required changed to move from
// the non-extensible factory to the exensible factory.  All other changes
// relative to the "factory" example are there to demonstrate new features.
/////////////////////////////////////////////////////////////////////////////
//non-extensible:  #include "G4PhysListFactory.hh"
#include "G4PhysListFactoryAlt.hh"
//use this for drop-in replacement:  using namespace g4alt;

/////////////////////////////////////////////////////////////////////////////
// headers needed to demonstrate new featues
/////////////////////////////////////////////////////////////////////////////

// allow ourselves to extend the short names for physics ctor addition/replace
// along the same lines as EMX, EMY, etc
#include "G4PhysListRegistry.hh"

// allow ourselves to give the user extra info about available physics ctors
#include "G4PhysicsConstructorFactory.hh"

// pull in a user defined physics list definition into the main program
// and register it with the factory (doesn't have to be the main program
// but the .o containing the declaration _must_ get linked/loaded)
#include "G4PhysListStamper.hh"  // defines macro for factory registration
#include "MySpecialPhysList.hh"
G4_DECLARE_PHYSLIST_FACTORY(MySpecialPhysList);

/////////////////////////////////////////////////////////////////////////////

#include "G4VModularPhysicsList.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {

  void PrintAvailable(G4int verbosity) {
    G4cout << G4endl;
    G4cout << "extensibleFactory: here are the available physics lists:"
           << G4endl;
    g4alt::G4PhysListFactory factory;
    factory.PrintAvailablePhysLists();

    // if user asked for extra verbosity then print physics ctors as well
    if ( verbosity > 1 ) {
      G4cout << G4endl;
      G4cout << "extensibleFactory: "
             << "here are the available physics ctors that can be added:"
             << G4endl;
      G4PhysicsConstructorRegistry* g4pctorFactory =
        G4PhysicsConstructorRegistry::Instance();
      g4pctorFactory->PrintAvailablePhysicsConstructors();
    }
  }

  void PrintUsage(G4int verbosity) {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " extensibleFactory [-m macro ] [-p physList ]"
           << " [-u UIsession] [-t nThreads]" << G4endl
           << " [-v | --verbose] [-h | --help]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
    G4cerr << "   note: -v can be repeated to increase verbosity." << G4endl;
    G4cerr << G4endl;

    if (verbosity>0) PrintAvailable(verbosity);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  if ( argc > 13 ) {
    PrintUsage(0);
    return 1;
  }

  G4String macro;
  G4String session;
  G4String physListName;
  char*    physListNameEnv = 0;
  G4String gdmlFileName;
#ifdef G4MULTITHREADED
  G4int nofThreads = 0;
#endif
  G4int verbosity = 0;

  for ( G4int i=1; i<argc; i=i+2 ) {
    G4String g4argv(argv[i]);  // convert only once
    if      ( g4argv == "-m" ) macro = argv[i+1];
    else if ( g4argv == "-u" ) session = argv[i+1];
    else if ( g4argv == "-p" ) physListName = argv[i+1];
#ifdef G4MULTITHREADED
    else if ( g4argv == "-t" ) {
      nofThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else if ( g4argv == "-v" || g4argv == "--verbose" ) {
      ++verbosity;  // verbose flag doesn't take an argument
      --i ;         // don't increment argc by two, just the one
    }
    else if ( g4argv == "-h" || g4argv == "--help" ) {
      PrintUsage(verbosity+1);
      return 1;
    }
    else {
      PrintUsage(0);
      return 1;
    }
  }

  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }

  // Choose the Random engine  //choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine());

  // Construct the run manager
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager();
  if ( nofThreads > 0 ) {
    runManager->SetNumberOfThreads(nofThreads);
  }
#else
  G4RunManager * runManager = new G4RunManager();
#endif

  // g4alt::G4PhysListFactoryAlt is the extensible factory
  // including the G4PhysListFactoryAlt.hh header and the line:
  //    using namespace g4alt;
  // would make this a drop-in replacement, but we'll list the explicit
  // namespace here just for clarity
  g4alt::G4PhysListFactory factory;
  G4VModularPhysicsList* physList = nullptr;

  // Show how an alternative default list could be set
  // but set the default to the normal default FTFP_BERT.
  // This is what is used when no -p flag is given and $PHYSLIST
  // is not defined in the environment.
  G4String defaultPhysListName = "FTFP_BERT";
  if ( verbosity > 0 ) {
    G4cout << "extensibleFactory: SetDefaultReferencePhysList to '"
           << defaultPhysListName << "' ('' = system default)"
           << G4endl << G4endl;
  }
  factory.SetDefaultReferencePhysList(defaultPhysListName);

  // set a short name for G4RadioactiveDecayPhysics
  G4PhysListRegistry* plreg = G4PhysListRegistry::Instance();
  plreg->AddPhysicsExtension("RADIO","G4RadioactiveDecayPhysics");
  plreg->AddPhysicsExtension("MYPHYSICS","MyG4PhysicsPhysics");
  if ( verbosity > 0 ) {
    G4cout << "extensibleFactory: adding extensions" << G4endl
           << "   RADIO     ===> G4RadioactiveDecayPhysics" << G4endl
           << "   MYPHYSICS ===> MyG4PhysicsPhysics" << G4endl
           << G4endl;
  }

  // Get Reference PhysicsList via its name, or if none given
  //    from environment varialb e$PHYSLIST, with fall back to a default
  if ( physListName.size() ) {
    if ( verbosity > 0 ) {
      G4cout << "extensibleFactory: explicitly using '"
             << physListName << "'" << G4endl;
    }
    physList = factory.GetReferencePhysList(physListName);
  } else {
    if ( verbosity > 0 ) {
      G4cout << "extensibleFactory: no -p flag;"
             << " using ReferencePhysList() ($PHYSLIST or default)" << G4endl;
    }
    physList = factory.ReferencePhysList();

    if ( ! physList ) {
      // failed?  get what the user set, but we couldn't find
      physListNameEnv = std::getenv("PHYSLIST");
      if ( physListNameEnv ) {
        G4cout << "extensibleFactory: $PHYSLIST="
               << physListNameEnv << G4endl;
      }
    }

  }

  // deal with failure to get what the user wanted
  // print what they _could_ use
  if ( ! physList ) {
    G4cerr << "extensibleFactory: PhysicsList '"
           << ( physListNameEnv ? physListNameEnv : physListName )
           << "' was not available in g4alt::PhysListFactory." << G4endl;
    PrintAvailable(verbosity);

    // if we can't get what the user asked for...
    //    don't go on to use something else, that's confusing
    G4ExceptionDescription ED;
    ED << "The factory for the physicslist ["
       << ( physListNameEnv ? physListNameEnv : physListName )
       << "] does not exist!"
       << G4endl;
    G4Exception("extensibleFactory",
                "extensibleFactory001", FatalException, ED);
    exit(42);
  }

  // Set mandatory initialization classes
  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(physList);

  // set user action classes
  ActionInitialization* actinit = 
    new ActionInitialization("extensibleFactory");
  runManager->SetUserInitialization(actinit);

// Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( macro.size() ) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else {
    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
