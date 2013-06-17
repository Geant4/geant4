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
// #ifdef G4MULTITHREADED
// #include "G4MTRunManager.hh"
// #else
#include "G4RunManager.hh"
// #endif
 
#include "G4UImanager.hh" 
#include "Tst68DetectorConstruction.hh"
#include "Tst68ActionInitialization.hh"
#include "FTFP_BERT.hh"
#include "Randomize.hh" 
#include <ctime>


int main(int argc,char** argv) { 

  CLHEP::Ranlux64Engine defaultEngine( 1234567, 4 ); 
 
  G4Random::setTheEngine( &defaultEngine ); 
  G4int seed = std::time( NULL ); 
  G4Random::setTheSeed( seed ); 
  G4cout << G4endl 
         << " ===================================================== " << G4endl 
         << " Initial seed = " << seed << G4endl 
	 << " ===================================================== " << G4endl 
	 << G4endl;

// #ifdef G4MULTITHREADED
//   G4MTRunManager* runManager = new G4MTRunManager;
//   runManager->SetNumberOfThreads(1);
// #else  
  G4RunManager* runManager = new G4RunManager;
// #endif
 
  runManager->SetUserInitialization( new Tst68DetectorConstruction ); 

  FTFP_BERT *thePL = new FTFP_BERT;

  //thePL->SetDefaultCutValue( 0.020 *mm ); // 20 microns 
  runManager->SetUserInitialization( thePL );
  runManager->SetUserInitialization(new Tst68ActionInitialization);
  runManager->Initialize(); 

  G4UImanager* UI = G4UImanager::GetUIpointer(); 
  if ( argc==1 ) {   // Define UI session for interactive mode. 
  } else {   // Batch mode 
    G4String command = "/control/execute "; 
    G4String fileName = argv[1]; 
    UI->ApplyCommand(command+fileName); 
  } 

  // job termination 
  delete runManager; 
  return 0; 
} 
