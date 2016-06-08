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
// $Id: exampleB04.cc,v 1.5 2002/05/31 11:46:23 dressel Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleB04
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"

#include "B04DetectorConstruction.hh"
#include "B04PhysicsList.hh"
#include "B04PrimaryGeneratorAction.hh"

#include "B04ImportanceDetectorConstruction.hh"

// Files specific for biasing
#include "G4ParallelImportanceSampler.hh"

int main(int argc, char **argv)
{  
  G4int numberOfEvent = 1000;

  G4String random_status_out_file, random_status_in_file;
  G4long myseed = 345354;

  HepRandom::setTheSeed(myseed);

  G4RunManager *runManager = new G4RunManager;
  
  // create the   "tracking"    detector      -----------------
  runManager->SetUserInitialization(new B04DetectorConstruction);
  //  ---------------------------------------------------
  runManager->SetUserInitialization(new B04PhysicsList);
  runManager->SetUserAction(new B04PrimaryGeneratorAction);
  runManager->Initialize();

  // create the detector for the importances
  B04ImportanceDetectorConstruction importancedetector;
  importancedetector.Construct();
  // the IStore is filled during detector construction
  G4VIStore &aIstore = *importancedetector.GetIStore();
  
  // create importance sampler to importance sample neutrons
  // in the importance geometry
  G4ParallelImportanceSampler pmgr(aIstore, "neutron");
  pmgr.Initialize();
  
  runManager->BeamOn(numberOfEvent);

  return 0;
}
