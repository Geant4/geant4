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
// $Id: G4PImportanceWWindowScoreSampler.cc,v 1.1 2002-05-31 10:16:02 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4PImportanceWWindowScoreSampler.cc
//
// ----------------------------------------------------------------------

#include "G4ParallelWorld.hh"
#include "G4PImportanceWWindowScoreSampler.hh"
#include "G4ParallelImportanceSampler.hh"
#include "G4ParallelWorld.hh"
#include "G4ProcessPlacer.hh"
#include "G4VIStore.hh"
#include "G4VPScorer.hh"
#include "G4PScoreProcess.hh"
#include "G4ParallelWeightWindowProcess.hh"
#include "G4WeightWindowAlgorithm.hh"

G4PImportanceWWindowScoreSampler::
G4PImportanceWWindowScoreSampler(G4VIStore &is, 
				 G4VPScorer &ascorer,
				 const G4String &particlename,
				 G4VWeightWindowAlgorithm &wwalg,
				 const G4VImportanceAlgorithm *ialg) : 
  fParticleName(particlename),
  fParallelWorld(*(new G4ParallelWorld(is.GetWorldVolume()))),
  fParallelImportanceSampler(*(new 
			       G4ParallelImportanceSampler(is,
							   fParticleName,
							   fParallelWorld,
							   ialg))),
   fPScorer(ascorer),
   fIstore(is),
   fPScoreProcess(0),
   fPWeightWindowProcess(0),
   fWWAlgorithm(wwalg)
{}

G4PImportanceWWindowScoreSampler::~G4PImportanceWWindowScoreSampler()
{
  if (fPScoreProcess) {
    G4ProcessPlacer placer(fParticleName);
    placer.RemoveProcess(fPScoreProcess);
    delete  fPScoreProcess;
  }
  if (fPWeightWindowProcess) {
    G4ProcessPlacer placer(fParticleName);
    placer.RemoveProcess(fPWeightWindowProcess);
    delete fPWeightWindowProcess;
  }
  delete &fParallelImportanceSampler;
  delete &fParallelWorld;
}

G4PScoreProcess *
G4PImportanceWWindowScoreSampler::CreateParallelScoreProcess()
{
  if (!fPScoreProcess) {
    fPScoreProcess = 
      new G4PScoreProcess(fParallelWorld.GetParallelStepper(), 
			  fPScorer);
  }
  return fPScoreProcess;
}

G4VProcess *
G4PImportanceWWindowScoreSampler::CreateWeightWindowProcess()
{
  if (!fPWeightWindowProcess) {
    fPWeightWindowProcess = 
      new G4ParallelWeightWindowProcess(fIstore,
					fParallelWorld.
					GetParallelStepper(),
					fWWAlgorithm);
  }
  return fPWeightWindowProcess;
}

void G4PImportanceWWindowScoreSampler::Initialize()
{
  G4ProcessPlacer placer(fParticleName);
  placer.AddProcessAsSecondDoIt(CreateParallelScoreProcess());
  fParallelImportanceSampler.Initialize();
  placer.AddProcessAsSecondDoIt(CreateWeightWindowProcess());
}
