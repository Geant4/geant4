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
// $Id: G4ParallelImportanceManager.cc,v 1.6 2002-05-30 12:55:42 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelImportanceManager.cc
//
// ----------------------------------------------------------------------

#include "G4ParallelImportanceManager.hh"
#include "G4ParallelManager.hh"
#include "G4ParallelWorld.hh"
#include "G4IStore.hh"
#include "G4ImportanceSampler.hh"
#include "G4ParallelImportanceProcess.hh"
#include "G4ProcessPlacer.hh"
#include "G4ImportanceAlgorithm.hh"

G4ParallelImportanceManager::
G4ParallelImportanceManager(G4VIStore &is,
			    const G4String &particlename,
			    const G4VImportanceAlgorithm *ialg)
 : fParallelManager(*(new G4ParallelManager(is.GetWorldVolume(), particlename))),
   fCreatedPM(true),
   fDeleteAlg( ( ! ialg) ),
   fIalgorithm(( (fDeleteAlg) ? 
		 new G4ImportanceAlgorithm : ialg)),
   fSampler(new G4ImportanceSampler(*fIalgorithm, 
                                    fParallelManager.
                                    GetParallelWorld().GetParallelStepper(),  
                                    is)),
   fParallelImportanceProcess(0)
{}

G4ParallelImportanceManager::
G4ParallelImportanceManager(G4VIStore &is, 
			    G4ParallelManager &pmanager,
			    const G4VImportanceAlgorithm *ialg)
 : fParallelManager(pmanager),
   fCreatedPM(false),
   fDeleteAlg( ( ! ialg) ),
   fIalgorithm(( (fDeleteAlg) ? 
		 new G4ImportanceAlgorithm : ialg)),
   fSampler(new G4ImportanceSampler(*fIalgorithm, 
                                    fParallelManager.
                                    GetParallelWorld().GetParallelStepper(),  
                                    is)),
   fParallelImportanceProcess(0)
{}
  
  
  
G4ParallelImportanceProcess *
G4ParallelImportanceManager::CreateParallelImportanceProcess()
{
  if (!fParallelImportanceProcess) {
    fParallelImportanceProcess = 
      new G4ParallelImportanceProcess(*fSampler, 
				      fParallelManager.
				      GetParallelWorld().
				      GetGeoDriver(), 
				      fParallelManager.
				      GetParallelWorld().
				      GetParallelStepper());
  }
  return fParallelImportanceProcess;
}

void G4ParallelImportanceManager::Initialize()
{
  G4ProcessPlacer placer(fParallelManager.GetParticleName());
  placer.AddProcessAsSecondDoIt(CreateParallelImportanceProcess());
}

G4ParallelImportanceManager::~G4ParallelImportanceManager()
{
  if (fParallelImportanceProcess) {
    G4ProcessPlacer placer(fParallelManager.GetParticleName());
    placer.RemoveProcess(fParallelImportanceProcess);
    delete fParallelImportanceProcess;
  }
  if (fCreatedPM) delete &fParallelManager;
  if (fDeleteAlg) delete fIalgorithm;
  delete fSampler;
}
