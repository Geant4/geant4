//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: exampleB05.cc,v 1.6 2002/05/31 11:46:23 dressel Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleB05
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "g4std/iostream"

#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "B05DetectorConstruction.hh"
#include "B05PhysicsList.hh"
#include "B05PrimaryGeneratorAction.hh"

// Files specific for biasing and scoring
#include "G4PIScorer.hh"
#include "G4Sigma.hh"
#include "G4MassImportanceScoreSampler.hh"

int main(int argc, char **argv)
{  
  G4std::ostream *myout = &G4cout;
  G4int numberOfEvent = 1000;

  G4String random_status_out_file, random_status_in_file;
  G4long myseed = 345354;

  HepRandom::setTheSeed(myseed);

  G4RunManager *runManager = new G4RunManager;
  
  // create the detector      ---------------------------
  B05DetectorConstruction *detector = new B05DetectorConstruction();
  runManager->SetUserInitialization(detector);
  //  ---------------------------------------------------
  runManager->SetUserInitialization(new B05PhysicsList);
  runManager->SetUserAction(new B05PrimaryGeneratorAction);
  runManager->Initialize();

  // the IStore is filled during detector construction
  G4VIStore &aIstore = *detector->GetIStore();

  // create the importance and scoring sampler for biasing and scoring 
  // in the tracking world
  G4PIScorer aScorer(aIstore); 
  G4MassImportanceScoreSampler mim(aIstore, aScorer, "neutron");
  mim.Initialize();

  G4UImanager* UI;

  UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/control/execute init.mac");   

  runManager->BeamOn(numberOfEvent);

  // print all the numbers calculated from the scorer
  *myout << "output aScorer, tracking geometry, neutron" << G4endl;
  *myout << aScorer << G4endl;
  *myout << "----------------------------------------------"  << G4endl;

  return 0;
}
