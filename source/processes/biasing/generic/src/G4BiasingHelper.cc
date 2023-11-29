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
#include "G4BiasingHelper.hh"

#include "G4ProcessManager.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4ParallelGeometriesLimiterProcess.hh"

G4bool G4BiasingHelper::ActivatePhysicsBiasing(G4ProcessManager* pmanager,
					       G4String physicsProcessToBias,
					       G4String wrappedName)
{
  G4VProcess* physicsProcess(0);
  
  G4ProcessVector* vprocess = pmanager->GetProcessList();
  for (G4int ip = 0 ; ip < (G4int)vprocess->size() ; ++ip)
    {
      if ( (*vprocess)[ip]->GetProcessName() == physicsProcessToBias )
	{
	  physicsProcess = (*vprocess)[ip];
	  break;
	}
    }
  
  // -- process not found, return "false" to tell about failure
  if ( physicsProcess == 0 ) return false;
  
  // -- process is not a physics one, return "false" to tell about failure
  G4int processType = physicsProcess->GetProcessType();
  if ( ( processType != 2 ) &&  // EM
       ( processType != 3 ) &&  // Optical
       ( processType != 4 ) &&  // Hadronic
       ( processType != 6 ) )   // Decay
    return false;

  // -- prevent wrapper of wrapper...
  if ( dynamic_cast< G4BiasingProcessInterface* >( physicsProcess ) ) return false;

  // -- remember process indeces:
  G4int    atRestIndex = pmanager->GetProcessOrdering(physicsProcess, idxAtRest   );
  G4int alongStepIndex = pmanager->GetProcessOrdering(physicsProcess, idxAlongStep);
  G4int  postStepIndex = pmanager->GetProcessOrdering(physicsProcess, idxPostStep );

  // -- now remove the physic process, that will be replaced by a wrapped version:
  G4VProcess* removed = pmanager->RemoveProcess(physicsProcess);
  if ( removed != physicsProcess )
    {
      G4ExceptionDescription ed;
      ed << "Internal inconsistency in processes handling. Please report !" << G4endl;
      G4Exception("G4BiasingHelper::ActivatePhysicsBiasing(...)",
		  "BIAS.GEN.01",
		  FatalException,
		  ed);
    }
  
  G4BiasingProcessInterface* biasingWrapper = new G4BiasingProcessInterface( physicsProcess,
									     atRestIndex    != ordInActive,
									     alongStepIndex != ordInActive,
									     postStepIndex  != ordInActive,
									     wrappedName);

  if ( alongStepIndex == -1 ) alongStepIndex = ordDefault;
  
  pmanager->AddProcess( biasingWrapper,
			atRestIndex,
			alongStepIndex,
			postStepIndex);

  return true;
}

void G4BiasingHelper::ActivateNonPhysicsBiasing(G4ProcessManager* pmanager,
						G4String nonPhysicsProcessName )
{
  G4BiasingProcessInterface* biasingNonPhys(nullptr);
  if ( nonPhysicsProcessName == "" ) biasingNonPhys = new G4BiasingProcessInterface();
  else                               biasingNonPhys = new G4BiasingProcessInterface(nonPhysicsProcessName );
  pmanager->AddProcess( biasingNonPhys,
			ordInActive,
			ordInActive,
			ordDefault);
}

G4ParallelGeometriesLimiterProcess* G4BiasingHelper::AddLimiterProcess(G4ProcessManager* pmanager, const G4String& processName)
{
  G4ParallelGeometriesLimiterProcess* toReturn = nullptr;
  
  G4ProcessVector* processList = pmanager->GetProcessList();
  G4bool noInstance = true;
  for (G4int i = 0 ; i < (G4int)processList->size() ; ++i)
    {
      G4VProcess* process = (*processList)[i];
      if ( dynamic_cast< G4ParallelGeometriesLimiterProcess* >( process ) )
	{
	  noInstance = false;
	  
	  G4ExceptionDescription ed;
	  ed << "Trying to re-add a G4ParallelGeometriesLimiterProcess process to the process manager for '"<<
	    pmanager->GetParticleType()->GetParticleName() << " (PDG : " <<  pmanager->GetParticleType()->GetPDGEncoding() << " )"
	     << " while one is already present." << G4endl;
	  G4Exception("G4BiasingHelper::AddBiasingProcessLimiter(G4ProcessManager* pmanager)",
		      "BIAS.GEN.28",
		      JustWarning, ed,
		      "Call ignored.");
	  break;
	}
    }
  
  if ( noInstance )
    {
      G4ParallelGeometriesLimiterProcess* biasingLimiter = new G4ParallelGeometriesLimiterProcess(processName);
      pmanager->AddProcess                ( biasingLimiter );
      pmanager->SetProcessOrderingToSecond( biasingLimiter, idxAlongStep );
      pmanager->SetProcessOrderingToLast  ( biasingLimiter, idxPostStep  );
      
      toReturn = biasingLimiter;
    }
  
  return toReturn;
}
