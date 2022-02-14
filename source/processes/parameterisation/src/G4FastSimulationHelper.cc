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
#include "G4FastSimulationHelper.hh"

#include "G4ProcessManager.hh"
#include "G4FastSimulationManagerProcess.hh"

void G4FastSimulationHelper::ActivateFastSimulation(G4ProcessManager* pmanager, G4String parallelGeometryName )
{
  G4FastSimulationManagerProcess* fastSimProcess;
  if ( parallelGeometryName.empty() ) {
    fastSimProcess = new G4FastSimulationManagerProcess("fastSimProcess_massGeom");
    // -- For the parametrisation envelope belonging to the mass geometry case, the G4FastSimulationManagerProcess
    // -- is a PostStep process, and ordering does not matter:
    pmanager-> AddDiscreteProcess(fastSimProcess);
  }
  else {
    fastSimProcess = new G4FastSimulationManagerProcess("fastSimProcess_parallelGeom",parallelGeometryName);
  // -- For the parallel geometry case, the G4FastSimulationManagerProcessz
  // -- is an Along+PostStep process, and ordering matters:
  pmanager->AddProcess(fastSimProcess);
  pmanager->SetProcessOrdering(fastSimProcess, idxAlongStep, 1);
  }
  // If the parallel world
  // exists (with parallel world physics), e.g. for the sensitive detector.
  // In that case make sure fast simulation is the first process to be checked by the steppping manager
  // (highest ordering) so that user can kill the particle and/or deposit energy, ignoring other processes.
  // Otherwise the parallel world physics (which is a StronglyFroced process) will invoke a PostStepDoIt
  // on the same step, leading to e.g. duplicated energy deposits.

  // Register as the process with highest ordering so it is checked as the first one,
  // and since it is exclusively forced no other process will be considered (to be invoked).
  pmanager->SetProcessOrderingToLast(fastSimProcess, idxPostStep);
}
