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
/// \file Hadr10.cc
/// \brief Main program of the hadronic/Hadr10 example
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Threading.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh" 
#include "G4PhysListFactory.hh"
#include "DetectorConstruction.hh" 
#include "ActionInitialization.hh"
#include "CLHEP/Random/Ranlux64Engine.h" 
#include "CLHEP/Random/MTwistEngine.h"


int main( int argc, char** argv ) { 

  CLHEP::Ranlux64Engine defaultEngine( 1234567, 4 );  
  CLHEP::HepRandom::setTheEngine( &defaultEngine ); 
  G4int seed = time( NULL ); 
  CLHEP::HepRandom::setTheSeed( seed ); 
  G4cout << G4endl 
         << " ===================================================== " << G4endl 
         << " Initial seed = " << seed << G4endl 
         << " ===================================================== " << G4endl 
         << G4endl;

  auto* runManager = G4RunManagerFactory::CreateRunManager();
  //runManager->SetNumberOfThreads(4);
  
  DetectorConstruction* pDetectorInstance = new DetectorConstruction;  
  runManager->SetUserInitialization( pDetectorInstance ); 

  // Physics list factory: use the PHYSLIST environmental variable.
  G4PhysListFactory factory;
  G4VModularPhysicsList* thePL = factory.ReferencePhysList();   

  runManager->SetUserInitialization( thePL );
  runManager->SetUserInitialization( new ActionInitialization );

  G4UImanager* UI = G4UImanager::GetUIpointer(); 
  if ( argc == 1 ) {
    // Define UI session for interactive mode. 
  } else {
    // Batch mode 
    G4String command = "/control/execute "; 
    G4String fileName = argv[1]; 
    UI->ApplyCommand(command+fileName); 
  } 

  // job termination 
  delete runManager; 
  return 0; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
