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
// $Id: G4ParallelImportanceScoreManager.cc,v 1.5 2002-05-30 11:14:39 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelImportanceScoreManager.cc
//
// ----------------------------------------------------------------------

#include "G4ParallelManager.hh"
#include "G4ParallelImportanceScoreManager.hh"
#include "G4ParallelImportanceManager.hh"
#include "G4ParallelWorld.hh"
#include "G4ProcessPlacer.hh"
#include "G4VIStore.hh"
#include "G4VPScorer.hh"
#include "G4PScoreProcess.hh"

G4ParallelImportanceScoreManager::
G4ParallelImportanceScoreManager(G4VIStore &is, 
				 G4VPScorer &ascorer,
				 const G4String &particlename,
				 const G4VImportanceAlgorithm *ialg)
 : fParallelManager(*(new G4ParallelManager(is.GetWorldVolume(), particlename))),
   fParallelImportanceManager(*(new 
			       G4ParallelImportanceManager(is,
							   fParallelManager,
							   ialg))),
   fPScorer(ascorer),
   fPScoreProcess(0)
{}

G4ParallelImportanceScoreManager::~G4ParallelImportanceScoreManager()
{
  if (fPScoreProcess) {
    G4ProcessPlacer placer(fParallelManager.GetParticleName());
    placer.RemoveProcess(fPScoreProcess);
    delete  fPScoreProcess;
  }
  delete &fParallelImportanceManager;
  delete &fParallelManager;
}

G4PScoreProcess *
G4ParallelImportanceScoreManager::CreateParallelScoreProcess()
{
  if (!fPScoreProcess) {
    fPScoreProcess = 
      new G4PScoreProcess(fParallelManager.GetParallelWorld().
			  GetParallelStepper(), fPScorer);
  }
  return fPScoreProcess;
}

void G4ParallelImportanceScoreManager::Initialize()
{
  G4ProcessPlacer placer(fParallelManager.GetParticleName());
  placer.AddProcessAsSecondDoIt(CreateParallelScoreProcess());
  fParallelImportanceManager.Initialize();
}
