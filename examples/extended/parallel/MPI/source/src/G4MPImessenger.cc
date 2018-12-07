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
/// @file G4MPImessenger.cc
/// @brief Define MPI commands

#include "mpi.h"
#include <sstream>
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UImanager.hh"
#include "G4UIparameter.hh"
#include "G4MPImanager.hh"
#include "G4MPImessenger.hh"
#include "G4VMPIseedGenerator.hh"

// --------------------------------------------------------------------------
G4MPImessenger::G4MPImessenger()
  : G4UImessenger()
{
  // /mpi
  dir_ = new G4UIdirectory("/mpi/");
  dir_-> SetGuidance("MPI control commands");

  // /mpi/verbose
  verbose_ = new G4UIcmdWithAnInteger("/mpi/verbose", this);
  verbose_-> SetGuidance("Set verbose level.");
  verbose_-> SetParameterName("verbose", false, false);
  verbose_->SetRange("verbose>=0 && verbose<=1");

  // /mpi/status
  status_ = new G4UIcmdWithoutParameter("/mpi/status", this);
  status_-> SetGuidance( "Show mpi status.");

  // /mpi/execute
  execute_ = new G4UIcmdWithAString("/mpi/execute", this);
  execute_-> SetGuidance("Execute a macro file. (=/control/execute)");
  execute_-> SetParameterName("fileName", false, false);

  // /mpi/beamOn
  beam_on_ = new G4UIcommand("/mpi/beamOn", this);
  beam_on_-> SetGuidance("Start a parallel run w/ thread.");

  G4UIparameter* p1= new G4UIparameter("numberOfEvent", 'i', true);
  p1-> SetDefaultValue(1);
  p1-> SetParameterRange("numberOfEvent>=0");
  beam_on_-> SetParameter(p1);

  G4UIparameter* p2= new G4UIparameter("divide", 'b', true);
  p2-> SetDefaultValue(true);
  beam_on_-> SetParameter(p2);

  // /mpi/.beamOn
  dot_beam_on_ = new G4UIcommand("/mpi/.beamOn", this);
  dot_beam_on_-> SetGuidance("Start a parallel run w/o thread.");

  p1= new G4UIparameter("numberOfEvent", 'i', true);
  p1-> SetDefaultValue(1);
  p1-> SetParameterRange("numberOfEvent>=0");
  dot_beam_on_-> SetParameter(p1);

  p2= new G4UIparameter("divide", 'b', true);
  p2-> SetDefaultValue(true);
  dot_beam_on_-> SetParameter(p2);

  // /mpi/masterWeight
  master_weight_ = new G4UIcmdWithADouble("/mpi/masterWeight", this);
  master_weight_-> SetGuidance("Set weight for master node.");
  master_weight_-> SetParameterName("weight", false, false);
  master_weight_-> SetRange("weight>=0. && weight<=1.");

  // /mpi/showSeeds
  show_seeds_ = new G4UIcmdWithoutParameter("/mpi/showSeeds", this);
  show_seeds_-> SetGuidance("Show seeds of MPI nodes.");

  // /mpi/setMasterSeed
  set_master_seed_ = new G4UIcmdWithAnInteger("/mpi/setMasterSeed", this);
  set_master_seed_-> SetGuidance("Set a master seed for the seed generator.");
  set_master_seed_-> SetParameterName("seed", false, false);

  // /mpi/setSeed
  set_seed_ = new G4UIcommand("/mpi/setSeed", this);
  set_seed_-> SetGuidance("Set a seed for a specified node.");

  p1 = new G4UIparameter("node", 'i', false);
  p1-> SetParameterRange("node>=0");
  set_seed_-> SetParameter(p1);

  p2 = new G4UIparameter("seed", 'i', false);
  set_seed_-> SetParameter(p2);
}

// --------------------------------------------------------------------------
G4MPImessenger::~G4MPImessenger()
{
  delete verbose_;
  delete status_;
  delete execute_;
  delete beam_on_;
  delete dot_beam_on_;
  delete master_weight_;
  delete show_seeds_;
  delete set_master_seed_;
  delete set_seed_;

  delete dir_;
}

// --------------------------------------------------------------------------
void G4MPImessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if ( command == verbose_ ) { // /mpi/verbose
    G4int lv = verbose_-> GetNewIntValue(newValue);
    g4mpi_-> SetVerbose(lv);

  } else if ( command == status_ ) { // /mpi/status
    g4mpi_-> ShowStatus();

  } else if ( command == execute_ ) { // /mpi/execute
    G4UImanager* UI = G4UImanager::GetUIpointer();
    g4mpi_-> ExecuteMacroFile(UI-> FindMacroPath(newValue));

  } else if ( command == beam_on_ ) { // /mpi/beamOn
    std::istringstream is(newValue);
    G4int nevent;
    G4bool qdivide;
    is >> nevent >> qdivide;
    g4mpi_-> BeamOn(nevent, qdivide);

  } else if ( command == dot_beam_on_ ) { // /mpi/.beamOn
    std::istringstream is(newValue);
    G4int nevent;
    G4bool qdivide;
    is >> nevent >> qdivide;
    g4mpi_-> BeamOn(nevent, qdivide);

  } else if ( command == master_weight_ ) { // /mpi/masterWeight
    G4double weight= master_weight_-> GetNewDoubleValue(newValue);
    g4mpi_-> SetMasterWeight(weight);

  } else if ( command == show_seeds_ ) { // /mpi/showSeeds
    g4mpi_-> ShowSeeds();

  } else if ( command == set_master_seed_ ) { // /mpi/setMasterSeed
    std::istringstream is(newValue);
    G4long seed;
    is >> seed;
    g4mpi_-> GetSeedGenerator()-> SetMasterSeed(seed);
    g4mpi_-> DistributeSeeds();

  } else if ( command == set_seed_ ) { // /mpi/setSeed
    std::istringstream is(newValue);
    G4int inode;
    G4long seed;
    is >> inode >> seed;
    g4mpi_-> SetSeed(inode, seed);
  }

  return;
}

// --------------------------------------------------------------------------
G4String G4MPImessenger::GetCurrentValue(G4UIcommand* command)
{
 G4String cv;

 if ( command == verbose_ ) {
   cv = verbose_-> ConvertToString(g4mpi_->GetVerbose());
 } else if ( command == master_weight_ ) {
   cv= master_weight_-> ConvertToString(g4mpi_->GetMasterWeight());
 }

 return cv;
}
