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
// $Id: exampleB03.cc,v 1.6 2002-05-31 11:46:23 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleB03
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "B03DetectorConstruction.hh"
#include "B03PhysicsList.hh"
#include "B03PrimaryGeneratorAction.hh"

// Files specific for biasing
#include "G4MassImportanceSampler.hh"

int main(int argc, char **argv)
{  

  G4int numberOfEvent = 1000;

  G4String random_status_out_file, random_status_in_file;
  G4long myseed = 345354;

  HepRandom::setTheSeed(myseed);

  G4RunManager *runManager = new G4RunManager;
  
  // create the detector      ---------------------------
  B03DetectorConstruction *detector = new B03DetectorConstruction();
  runManager->SetUserInitialization(detector);
  //  ---------------------------------------------------
  runManager->SetUserInitialization(new B03PhysicsList);
  runManager->SetUserAction(new B03PrimaryGeneratorAction);
  runManager->Initialize();

  // the IStore is filled during detector construction
  G4VIStore &aIstore = *detector->GetIStore();
  // create the importance sampler for biasing in the tracking world
  G4MassImportanceSampler mim(aIstore, "neutron");
  mim.Initialize();

  G4UImanager* UI;

  UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/control/execute init.mac");   

  runManager->BeamOn(numberOfEvent);

  return 0;
}
