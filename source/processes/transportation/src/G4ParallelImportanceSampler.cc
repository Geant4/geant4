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
// $Id: G4ParallelImportanceSampler.cc,v 1.1 2002-05-31 10:16:02 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelImportanceSampler.cc
//
// ----------------------------------------------------------------------

#include "G4ParallelImportanceSampler.hh"
#include "G4ParallelWorld.hh"
#include "G4ParallelWorld.hh"
#include "G4IStore.hh"
#include "G4ImportanceSplitExaminer.hh"
#include "G4ParallelImportanceProcess.hh"
#include "G4ProcessPlacer.hh"
#include "G4ImportanceAlgorithm.hh"

G4ParallelImportanceSampler::
G4ParallelImportanceSampler(G4VIStore &is,
			    const G4String &particlename,
			    const G4VImportanceAlgorithm *ialg)
  :fParticleName(particlename), 
  fParallelWorld(*(new G4ParallelWorld(is.GetWorldVolume()))),
  fCreatedPW(true),
  fDeleteAlg( ( ! ialg) ),
  fIalgorithm(( (fDeleteAlg) ? 
		new G4ImportanceAlgorithm : ialg)),
  fExaminer(new G4ImportanceSplitExaminer(*fIalgorithm, 
					  fParallelWorld.
					  GetParallelStepper(),  
					  is)),
   fParallelImportanceProcess(0)
{}

G4ParallelImportanceSampler::
G4ParallelImportanceSampler(G4VIStore &is,
			    const G4String &pname,
			    G4ParallelWorld &pworld,
			    const G4VImportanceAlgorithm *ialg) : 
  fParticleName(pname), 
  fParallelWorld(pworld),
  fCreatedPW(false),
  fDeleteAlg( ( ! ialg) ),
  fIalgorithm(( (fDeleteAlg) ? 
		new G4ImportanceAlgorithm : ialg)),
  fExaminer(new G4ImportanceSplitExaminer(*fIalgorithm, 
					  fParallelWorld.
					  GetParallelStepper(),  
					  is)),
   fParallelImportanceProcess(0)
{}
  
  
  
G4ParallelImportanceProcess *
G4ParallelImportanceSampler::CreateParallelImportanceProcess()
{
  if (!fParallelImportanceProcess) {
    fParallelImportanceProcess = 
      new G4ParallelImportanceProcess(*fExaminer, 
				      fParallelWorld.
				      GetGeoDriver(), 
				      fParallelWorld.
				      GetParallelStepper());
  }
  return fParallelImportanceProcess;
}

void G4ParallelImportanceSampler::Initialize()
{
  G4ProcessPlacer placer(fParticleName);
  placer.AddProcessAsSecondDoIt(CreateParallelImportanceProcess());
}

G4ParallelImportanceSampler::~G4ParallelImportanceSampler()
{
  if (fParallelImportanceProcess) {
    G4ProcessPlacer placer(fParticleName);
    placer.RemoveProcess(fParallelImportanceProcess);
    delete fParallelImportanceProcess;
  }
  if (fCreatedPW) delete &fParallelWorld;
  if (fDeleteAlg) delete fIalgorithm;
  delete fExaminer;
}
