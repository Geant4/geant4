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
// $Id: G4ParallelScoreManager.cc,v 1.4 2002-05-24 08:17:20 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelScoreManager.cc
//
// ----------------------------------------------------------------------

#include "G4ParallelScoreManager.hh"
#include "G4ParallelManager.hh"
#include "G4PScoreProcess.hh"
#include "G4ProcessPlacer.hh"
#include "G4ParallelWorld.hh"

G4ParallelScoreManager::
G4ParallelScoreManager(G4VPhysicalVolume &worldvolume,
		       const G4String &particlename,
		       G4VPScorer &scorer)
 : fParallelManager(*(new G4ParallelManager(worldvolume, particlename))),
   fPScorer(scorer),
   fPScorerProcess(0)
{}

G4ParallelScoreManager::~G4ParallelScoreManager()
{
  if (fPScorerProcess) {
    G4ProcessPlacer placer(fParallelManager.GetParticleName());
    placer.RemoveProcess(fPScorerProcess);
    delete  fPScorerProcess;
  }
  delete &fParallelManager;
}


G4PScoreProcess *G4ParallelScoreManager::CreateParallelScoreProcess()
{
  if (!fPScorerProcess) {
    fPScorerProcess = 
      new G4PScoreProcess(fParallelManager.GetParallelWorld().
			  GetParallelStepper(), fPScorer);
  }
  return fPScorerProcess;
}

void G4ParallelScoreManager::Initialize()
{
  G4ProcessPlacer placer(fParallelManager.GetParticleName());
  placer.AddProcessAsSecondDoIt(CreateParallelScoreProcess());
  fParallelManager.Initialize();
}
