#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include <unistd.h>
#define _GNU_SOURCE
#include <getopt.h>
#include "g4std/set"

#include <iomanip>

#include "B02DetectorConstruction.hh"
#include "B02PhysicsList.hh"
#include "B02PrimaryGeneratorAction.hh"

#include "B02ScoringDetectorConstruction.hh"

// Files specific for scoring 
#include "G4PScorer.hh"
#include "G4Sigma.hh"
#include "G4ParallelScoreManager.hh"

int main(int argc, char **argv) {
  

  ostream *myout = &G4cout;
  int numberOfEvent = 1000;

  G4String random_status_out_file, random_status_in_file;
  G4long myseed = 345354;

  HepRandom::setTheSeed(myseed);

  G4RunManager *runManager = new G4RunManager;
  
  // create the   "tracking"    detector      -----------------
  runManager->
    SetUserInitialization(new B02DetectorConstruction);
  //  ---------------------------------------------------
  runManager->SetUserInitialization(new B02PhysicsList);
  runManager->SetUserAction(new B02PrimaryGeneratorAction);
  runManager->Initialize();

  // create the detector for scoring
  B02ScoringDetectorConstruction scoringdetector;
  
  // create scorer and manager to score gammas in the "scoring" detector
  G4PScorer pScorer;
  G4ParallelScoreManager pmgr(*(scoringdetector.Construct()), 
			      "gamma", pScorer);
  pmgr.Initialize();

  runManager->BeamOn(numberOfEvent);

  // ======= after running ============================

  // print all the numbers calculated from the scorer
  *myout << "output pScorer, scoring detector, gamma" << G4endl;
  *myout << pScorer << G4endl;
  *myout << "----------------------------------------------"  << G4endl;

  return 0;
}

