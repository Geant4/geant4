#include "G4BiasingHelper.hh"

#include "G4ProcessManager.hh"
#include "G4BiasingProcessInterface.hh"

G4bool G4BiasingHelper::ActivatePhysicsBiasing(G4ProcessManager* pmanager,
					       G4String physicsProcessToBias,
					       G4String wrappedName)
{
  G4VProcess* physicsProcess(0);
  
  G4ProcessVector* vprocess = pmanager->GetProcessList();
  for (G4int ip = 0 ; ip < vprocess->size() ; ip++)
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
  G4BiasingProcessInterface* biasingNonPhys(0);
  if ( nonPhysicsProcessName == "" ) biasingNonPhys = new G4BiasingProcessInterface();
  else                               biasingNonPhys = new G4BiasingProcessInterface(nonPhysicsProcessName );
  pmanager->AddProcess( biasingNonPhys,
			ordInActive,
			ordInActive,
			ordDefault);
}
