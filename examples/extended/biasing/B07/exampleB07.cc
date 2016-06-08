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
// $Id: exampleB07.cc,v 1.5 2002/05/31 11:46:24 dressel Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleB07
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "g4std/iostream"

#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"

#include "B07DetectorConstruction.hh"
#include "B07PhysicsList.hh"
#include "B07PrimaryGeneratorAction.hh"
#include "B07ImportanceDetectorConstruction.hh"

// Files specific for biasing and scoring
#include "G4ParallelImportanceScoreSampler.hh"
#include "G4PIScorer.hh"
#include "G4Sigma.hh"

int main(int argc, char **argv)
{  
  G4std::ostream *myout = &G4cout;
  G4int numberOfEvent = 10000;

  G4String random_status_out_file, random_status_in_file;
  G4long myseed = 345354;

  HepRandom::setTheSeed(myseed);

  G4RunManager *runManager = new G4RunManager;
  
  // create the   "tracking"    detector      -----------------
  runManager->SetUserInitialization(new B07DetectorConstruction);
  //  ---------------------------------------------------
  runManager->SetUserInitialization(new B07PhysicsList);
  runManager->SetUserAction(new B07PrimaryGeneratorAction);
  runManager->Initialize();

  // create the detector for the importances and scoring
  B07ImportanceDetectorConstruction importancedetector;
  importancedetector.Construct();
  // the IStore is filled during detector construction
  G4VIStore &aIstore = *importancedetector.GetIStore();
  
  // ceate two scorers 
  G4PIScorer nScorer(aIstore); // scorer for neutrons
  G4PIScorer gScorer(aIstore); // scorer for gammas
  
  // create two importance sampler to importance sample neutrons
  // and gammas in the importance geometry
  G4ParallelImportanceScoreSampler npmgr(aIstore, nScorer, "neutron");
  npmgr.Initialize();
  G4ParallelImportanceScoreSampler gpmgr(aIstore, gScorer, "gamma");
  gpmgr.Initialize();
  
  runManager->BeamOn(numberOfEvent);

  // print all the numbers calculated from the scorer
  *myout << "output nScorer, importance geometry, neutron" << G4endl;
  *myout << nScorer << G4endl;
  *myout << "----------------------------------------------"  << G4endl;

    // print all the numbers calculated from the scorer
  *myout << "output gScorer, importance geometry, gamma" << G4endl;
  *myout << gScorer << G4endl;
  *myout << "----------------------------------------------"  << G4endl;

  return 0;
}
