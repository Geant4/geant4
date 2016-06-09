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
// $Id: BrachySimulation.cc
// GEANT4 tag $Name: geant4-08-02 $
//
// --------------------------------------------------------------
//                 GEANT 4 - Brachytherapy example
// --------------------------------------------------------------
//
// Code developed by: S.Guatelli, K. Moscicki
//
//
//    *******************************
//    *                             *
//    *    BrachySimualtion.cc      *
//    *                             *
//    *******************************
//
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "BrachyFactoryIr.hh"
#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "BrachyEventAction.hh"
#include "BrachyDetectorConstruction.hh"
#include "BrachyPhysicsList.hh"
#include "BrachyPhantomSD.hh"
#include "BrachyPrimaryGeneratorActionIr.hh"
#include "G4SDManager.hh"
#include "BrachyRunAction.hh"
#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4UImessenger.hh"
#include "BrachySimulation.hh"

BrachySimulation::BrachySimulation(G4int sd)
{ 
  pRunManager = 0;
  seed = sd;
}

BrachySimulation::~BrachySimulation()

{ }
void BrachySimulation::setSeed(G4int sd)
{
  seed = sd;
  CLHEP::HepRandom::setTheSeed(seed);
}

G4bool BrachySimulation::initialize(int ,char** )
{ 
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
 
  G4cout << "G4 initializing" << G4endl;
 
  pRunManager = new G4RunManager;

  G4String sensitiveDetectorName = "Phantom";

  BrachyDetectorConstruction  *pDetectorConstruction = new  BrachyDetectorConstruction(sensitiveDetectorName);

  pRunManager -> SetUserInitialization(pDetectorConstruction);
  pRunManager -> SetUserInitialization(new BrachyPhysicsList);
  
  // output environment variables:
#ifdef G4ANALYSIS_USE
  G4cout << G4endl << G4endl << G4endl 
	 << " User Environment " << G4endl
	 << " Using AIDA 3.2.1 analysis " << G4endl;
# else
  G4cout << G4endl << G4endl << G4endl 
	 << " User Environment " << G4endl
	 << " G4ANALYSIS_USE environment variable not set, NO ANALYSIS " 
	 << G4endl;
#endif

  BrachyEventAction *pEventAction = new BrachyEventAction();
  pRunManager -> SetUserAction(pEventAction );

  BrachyRunAction *pRunAction = new BrachyRunAction();
  pRunManager -> SetUserAction(pRunAction);

  //Initialize G4 kernel
  pRunManager -> Initialize();

  return true;
}
void BrachySimulation::executeMacro(std::string macroFileName)
{
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  // Batch mode

  std::string fileName = macroFileName;
  G4cout << fileName << " <---------- in batch -----------------" << G4endl;
  G4cout << macroFileName << " executed"<<G4endl;

  std::string command = "/control/execute ";

  if(!UI)
    G4cout << "FATAL ERROR: UI pointer does not exist" << G4endl;
  else
    UI -> ApplyCommand(command+fileName);    
}
std::string  BrachySimulation::getOutputFilename()
{
  return "brachytherapy.xml";
}

void BrachySimulation::finish()
{
  delete pRunManager;
}

//// This is all the application needs to run in parallel mode through DIANE
extern "C" 
DIANE::IG4Simulation* createG4Simulation(int seed) 
{return new BrachySimulation(seed); }
///////////
