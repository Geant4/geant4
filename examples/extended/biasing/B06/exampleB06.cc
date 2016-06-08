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
// $Id: exampleB06.cc,v 1.7 2002/05/31 11:46:24 dressel Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleB06
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "g4std/iostream"

#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"

#include "B06DetectorConstruction.hh"
#include "B06PhysicsList.hh"
#include "B06PrimaryGeneratorAction.hh"
#include "B06ImportanceDetectorConstruction.hh"

// Files specific for biasing and scoring
#include "G4ParallelImportanceScoreSampler.hh"
#include "B06Scorer.hh"
#include "G4Sigma.hh"

// class for special output
#include "B06ScorePrinter.hh"

int main(int argc, char **argv)
{  
  G4std::ostream *myout = &G4cout;
  G4int numberOfEvent = 1000;

  G4String random_status_out_file, random_status_in_file;
  G4long myseed = 345354;

  HepRandom::setTheSeed(myseed);

  G4RunManager *runManager = new G4RunManager;
  
  // create the   "tracking"    detector      -----------------
  runManager->SetUserInitialization(new B06DetectorConstruction);
  //  ---------------------------------------------------
  runManager->SetUserInitialization(new B06PhysicsList);
  runManager->SetUserAction(new B06PrimaryGeneratorAction);
  runManager->Initialize();

  // create the detector for the importances and scoring
  B06ImportanceDetectorConstruction importancedetector;
  // the IStore is filled during detector construction
  G4VIStore &aIstore = *importancedetector.GetIStore();
  
  // create importance sampler to importance sample neutrons
  // in the importance geometry
  B06Scorer iScorer(aIstore);
  G4ParallelImportanceScoreSampler pmgr(aIstore, iScorer, "neutron");
  pmgr.Initialize();
  
  runManager->BeamOn(numberOfEvent);

    // print all the numbers calculated from the scorer
  *myout << "output iScorer, importance geometry, neutron" << G4endl;
  *myout << iScorer << G4endl;
  *myout << "----------------------------------------------"  << G4endl;

  // special output

  B06ScorePrinter sp;
  sp.PrintHeader(myout);
  sp.PrintTable(iScorer.GetMapPtkTallys(), myout);

  return 0;
}
