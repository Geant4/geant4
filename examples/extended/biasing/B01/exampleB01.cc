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
// $Id: exampleB01.cc,v 1.12 2002-08-29 15:37:24 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleB01
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

#include "B01DetectorConstruction.hh"
#include "B01PhysicsList.hh"
#include "B01PrimaryGeneratorAction.hh"

// Files specific for scoring 
#include "G4Scorer.hh"
#include "G4MassScoreSampler.hh"

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
  
  // create the detector      ---------------------------
  runManager->SetUserInitialization(new B01DetectorConstruction);
  //  ---------------------------------------------------
  runManager->SetUserInitialization(new B01PhysicsList);
  runManager->SetUserAction(new B01PrimaryGeneratorAction);
  runManager->Initialize();

  // create scorer and sampler to score neutrons in the detector

  G4Scorer scorer; // a scorer 

  G4MassScoreSampler msm(scorer, "neutron"); // to be done after 
  msm.Initialize();                           // runManager->Initialize()

  runManager->BeamOn(numberOfEvent);

  // ======= after running ============================

  G4ScoreTable sp;
  sp.Print(scorer.GetMapGeometryCellCellScorer(), myout);

  return 0;
}
