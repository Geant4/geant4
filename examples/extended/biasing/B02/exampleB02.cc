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
// $Id: exampleB02.cc,v 1.11 2002-08-29 15:37:25 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleB02
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "g4std/set"
#include "g4std/iomanip"

#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"

#include "B02DetectorConstruction.hh"
#include "B02PhysicsList.hh"
#include "B02PrimaryGeneratorAction.hh"

#include "B02ScoringDetectorConstruction.hh"

// Files specific for scoring
#include "G4Scorer.hh"
#include "G4ParallelScoreSampler.hh"

// table for scores
#include "G4ScoreTable.hh"


int main(int argc, char **argv)
{  

  G4std::ostream *myout = &G4cout;
  G4int numberOfEvent = 1000;

  G4String random_status_out_file, random_status_in_file;
  G4long myseed = 345354;

  HepRandom::setTheSeed(myseed);

  G4RunManager *runManager = new G4RunManager;
  
  // create the   "tracking"    detector      -----------------
  runManager->SetUserInitialization(new B02DetectorConstruction);
  //  ---------------------------------------------------
  runManager->SetUserInitialization(new B02PhysicsList);
  runManager->SetUserAction(new B02PrimaryGeneratorAction);
  runManager->Initialize();

  // create the detector for scoring
  B02ScoringDetectorConstruction scoringdetector;
  
  // create scorer and sampler to score neutrons in the "scoring" detector

  G4Scorer scorer; // a scorer   

  G4ParallelScoreSampler pmgr(*(scoringdetector.Construct()), 
			      "neutron", scorer); // sample w. scoting
  pmgr.Initialize();      // set up the necessary processes

  runManager->BeamOn(numberOfEvent);

  // ======= after running ============================

  // print a table
  G4ScoreTable sp;
  sp.Print(scorer.GetMapGeometryCellCellScorer(), myout);

 

  return 0;
}

