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
// $Id: G4MPImessenger.cc,v 1.3 2010-12-03 08:21:29 kmura Exp $
// $Name: not supported by cvs2svn $
//
// ====================================================================
//   G4MPImessenger.cc
//
//                                         2007 Q
// ====================================================================
#include "G4MPImessenger.hh"
#include "G4MPImanager.hh"
#include "G4VMPIseedGenerator.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include <sstream>

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////////////////////////////////
G4MPImessenger::G4MPImessenger( G4MPImanager* manager)
  : G4UImessenger(), g4MPI(manager)
//////////////////////////////////////////////////////
{
  // /mpi
  dir= new G4UIdirectory("/mpi/");
  dir-> SetGuidance("MPI control commands");

  // /mpi/verbose
  verbose= new G4UIcmdWithAnInteger("/mpi/verbose", this);
  verbose-> SetGuidance("Set verbose level.");
  verbose-> SetParameterName("verbose", false, false);
  verbose->SetRange("verbose>=0 && verbose<=1");

  // /mpi/status
  status= new G4UIcmdWithoutParameter("/mpi/status", this);
  status-> SetGuidance( "Show mpi status.");

  // /mpi/execute
  execute= new G4UIcmdWithAString("/mpi/execute", this);
  execute-> SetGuidance("Execute a macro file. (=/control/execute)");
  execute-> SetParameterName("fileName", false, false);

  // /mpi/beamOn
  beamOn= new G4UIcommand("/mpi/beamOn", this);
  beamOn-> SetGuidance("Start a parallel run w/ thread.");

  G4UIparameter* p1= new G4UIparameter("numberOfEvent", 'i', true);
  p1-> SetDefaultValue(1);
  p1-> SetParameterRange("numberOfEvent>=0");
  beamOn-> SetParameter(p1);

  G4UIparameter* p2= new G4UIparameter("divide", 'b', true);
  p2-> SetDefaultValue(true);
  beamOn-> SetParameter(p2);

  // /mpi/.beamOn
  dotbeamOn= new G4UIcommand("/mpi/.beamOn", this);
  dotbeamOn-> SetGuidance("Start a parallel run w/o thread.");

  p1= new G4UIparameter("numberOfEvent", 'i', true);
  p1-> SetDefaultValue(1);
  p1-> SetParameterRange("numberOfEvent>=0");
  dotbeamOn-> SetParameter(p1);

  p2= new G4UIparameter("divide", 'b', true);
  p2-> SetDefaultValue(true);
  dotbeamOn-> SetParameter(p2);

  // /mpi/weightForMaster
  masterWeight= new G4UIcmdWithADouble("/mpi/masterWeight", this);
  masterWeight-> SetGuidance("Set weight for master node.");
  masterWeight-> SetParameterName("weight", false, false);
  masterWeight-> SetRange("weight>=0. && weight<=1.");

  // /mpi/showSeeds
  showSeeds= new G4UIcmdWithoutParameter("/mpi/showSeeds", this);
  showSeeds-> SetGuidance("Show seeds of MPI nodes.");

  // /mpi/setMasterSeed
  setMasterSeed= new G4UIcmdWithAnInteger("/mpi/setMasterSeed", this);
  setMasterSeed-> SetGuidance("Set a master seed for the seed generator.");
  setMasterSeed-> SetParameterName("seed", false, false);

  // /mpi/setSeed
  setSeed= new G4UIcommand("/mpi/setSeed", this);
  setSeed-> SetGuidance("Set a seed for a specified node.");

  p1= new G4UIparameter("node", 'i', false);
  p1-> SetParameterRange("node>=0");
  setSeed-> SetParameter(p1);

  p2= new G4UIparameter("seed", 'i', false);
  setSeed-> SetParameter(p2);
}


/////////////////////////////////
G4MPImessenger::~G4MPImessenger()
/////////////////////////////////
{
  delete verbose;
  delete status;
  delete execute;
  delete beamOn;
  delete dotbeamOn;
  delete masterWeight;
  delete showSeeds;
  delete setMasterSeed;
  delete setSeed;

  delete dir;
}


//////////////////////////////////////////////////////////////////////////
void G4MPImessenger::SetNewValue( G4UIcommand* command, G4String newValue)
//////////////////////////////////////////////////////////////////////////
{
  if (command == verbose) { // /mpi/verbose
    G4int lv= verbose-> GetNewIntValue(newValue);
    g4MPI-> SetVerbose(lv);

  } else if (command == status){ // /mpi/status
    g4MPI-> ShowStatus();

  } else if (command == execute){ // /mpi/execute
    g4MPI-> ExecuteMacroFile(newValue);

  } else if (command == beamOn){ // /mpi/beamOn
    std::istringstream is(newValue);
    G4int nevent;
    G4bool qdivide;
    is >> nevent >> qdivide;
    g4MPI-> BeamOn(nevent, qdivide);

  } else if (command == dotbeamOn){ // /mpi/.beamOn
    std::istringstream is(newValue);
    G4int nevent;
    G4bool qdivide;
    is >> nevent >> qdivide;
    g4MPI-> BeamOn(nevent, qdivide);

  } else if (command == masterWeight){ // /mpi/masterWeight
    G4double weight= masterWeight-> GetNewDoubleValue(newValue);
    g4MPI-> SetMasterWeight(weight);

  } else if (command == showSeeds){ // /mpi/showSeeds
    g4MPI-> ShowSeeds();

  } else if (command == setMasterSeed){ // /mpi/setMasterSeed
    std::istringstream is(newValue);
    G4long seed;
    is >> seed;
    g4MPI-> GetSeedGenerator()-> SetMasterSeed(seed);
    g4MPI-> DistributeSeeds();

  } else if (command == setSeed){ // /mpi/setSeed
    std::istringstream is(newValue);
    G4int inode;
    G4long seed;
    is >> inode >> seed;
    g4MPI-> SetSeed(inode, seed);
  }

  return;
}


/////////////////////////////////////////////////////////////
G4String G4MPImessenger::GetCurrentValue(G4UIcommand* command)
/////////////////////////////////////////////////////////////
{
 G4String cv;

 if (command == verbose) {
   cv= verbose-> ConvertToString(g4MPI->GetVerbose());
 } else if (command == masterWeight) {
   cv= masterWeight-> ConvertToString(g4MPI->GetMasterWeight());
 }

 return cv;
}

