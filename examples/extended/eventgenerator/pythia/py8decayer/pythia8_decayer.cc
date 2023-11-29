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
/// \file eventgenerator/pythia/py8decayer/pythia8_decayer.cc
/// \brief Main program of the pythia8_decayer example

#include "G4PhysListFactoryAlt.hh"
#include "G4PhysListRegistry.hh"

#include "Py8DecayerPhysics.hh"

#include "DetConstruction.hh"
#include "SingleParticleGun.hh"
#include "ActionInitialization.hh"

#include "G4RunManagerFactory.hh"
#include "FTFP_BERT.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4GeometryManager.hh"

#include "G4PhysListFactoryAlt.hh"

// allow ourselves to extend the short names for physics ctor addition/replace
// along the same lines as EMX, EMY, etc
#include "G4PhysListRegistry.hh"

// allow ourselves to give the user extra info about available physics ctors
#include "G4PhysicsConstructorFactory.hh"

// forward declaration
void PrintAvailable(G4int verb=1);

int main(int argc,char** argv)
{

   // pick physics list
   std::string physListName = "FTFP_BERT+PY8DK";

   for ( G4int i=1; i<argc; i=i+2 ) {
      G4String g4argv(argv[i]);  // convert only once
      if ( g4argv == "-p" ) physListName = argv[i+1];
   }

   // Choose the Random engine
   //
   G4Random::setTheEngine(new CLHEP::RanecuEngine);

   // Construct the default run manager
   //
   auto runManager =
      G4RunManagerFactory::CreateRunManager( G4RunManagerType::SerialOnly );

   // NOTE: in practice, it will only work in MT mode
   //
   runManager->SetNumberOfThreads(5);

   // Set mandatory initialization classes
   //
   // Geometry
   //
   runManager->SetUserInitialization(new DetConstruction);
   //
   // Physics
   //
   //G4VModularPhysicsList* physicsList = new FTFP_BERT();
   //physicsList->RegisterPhysics(new Py8DecayerPhysics());

   g4alt::G4PhysListFactory plFactory;
   G4VModularPhysicsList* physicsList = nullptr;
   plFactory.SetDefaultReferencePhysList("NO_DEFAULT_PHYSLIST");

   // set a short name for the plugin
   G4PhysListRegistry* plReg = G4PhysListRegistry::Instance();
   plReg->AddPhysicsExtension("PY8DK","Py8DecayerPhysics");

   physicsList = plFactory.GetReferencePhysList(physListName);

   if ( ! physicsList ) {
     PrintAvailable( 1 );

     // if we can't get what the user asked for...
     //    don't go on to use something else, that's confusing
     G4ExceptionDescription ed;
     ed << "The factory for the physicslist [" 
        << physListName << "] does not exist!" 
        << G4endl;
     G4Exception("extensibleFactory",
                 "extensibleFactory001", FatalException, ed);
     exit(42);
   }

   runManager->SetUserInitialization(physicsList);
   //
   // Set user action classes, e.g. prim.generator (tau- gun), etc.
   // NOTE: setting particle gun (prim. generator) is done in
   //       ActionInit class (see src/ActionInitialization.cc)
   //
   runManager->SetUserInitialization( new ActionInitialization() );

   // Run initialization
   //
   runManager->GeometryHasBeenModified();
   runManager->InitializeGeometry();
   runManager->InitializePhysics();
   runManager->Initialize();
   //
   runManager->ConfirmBeamOnCondition();
   runManager->RunInitialization(); 

   // Run event loop (explicitly rather than via BeamOn)
   //
   for (G4int iev=0; iev<5; ++iev)
   {
      runManager->ProcessOneEvent( iev );
   }

   // It speaks for itself ^_^
   //
   runManager->RunTermination();
    
   // Job termination
   // Free the store: user actions, physics_list and detector_description are
   //                 owned and deleted by the run manager, so they should not
   //                 be deleted in the main() program !
   delete runManager;

   return 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
