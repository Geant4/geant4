#include "GB01BOptrChangeCrossSection.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4BOptnChangeCrossSection.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4VProcess.hh"

#include "Randomize.hh"

GB01BOptrChangeCrossSection::GB01BOptrChangeCrossSection(G4String particleName,
                                                         G4String         name)
  : G4VBiasingOperator(name),
    fFirstProcess(0), 
    fLastProcess(0),
    fSetup(true)
{
  fParticleToBias = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
  // --> put a G4Exception
  if ( fParticleToBias == 0 ) G4cout << " ********* particle not found ****** " << G4endl;
  
}

GB01BOptrChangeCrossSection::~GB01BOptrChangeCrossSection()
{
  for ( std::map< const G4BiasingProcessInterface*, G4BOptnChangeCrossSection* >::iterator 
          it = fChangeCrossSectionOperations.begin() ;
        it != fChangeCrossSectionOperations.end() ;
        it++ ) delete (*it).second;
}

G4VBiasingOperation* 
GB01BOptrChangeCrossSection::ProposeOccurenceBiasingOperation(const G4Track*            track, 
                                                                 const G4BiasingProcessInterface*
                                                                                  callingProcess)
{

  // -- Check if current particle type is the one to bias:
  if ( track->GetDefinition() != fParticleToBias ) return 0;

  // ---------------
  // -- Setup stage:
  // ---------------
  // -- ( Might consider moving this in a less called method. )
  // -- Start by remembering processes under biasing, create needed biasing operations
  // -- and assocation operations to processes:
  if ( fSetup )
    {
      if ( ( fFirstProcess == 0 ) && 
           ( callingProcess->GetIsFirstPostStepGPILInterface() ) ) fFirstProcess = callingProcess;
      if ( fLastProcess == 0 )
        {
          G4String operationName = "XSchange-" + 
            callingProcess->GetWrappedProcess()->GetProcessName();
          fChangeCrossSectionOperations[callingProcess] = 
            new G4BOptnChangeCrossSection(operationName);
          if ( callingProcess->GetIsLastPostStepGPILInterface() )
            {
              fLastProcess = callingProcess;
              fSetup = false;
            }
        }
    }


  // -- Check if the analog cross-section undefined : this is for example the case of conversion
  // -- when a gamma is below e+e- creation threshold. Does nothing (ie, let analog
  // -- process to deal with the case)
  G4double analogInteractionLength =  
    callingProcess->GetWrappedProcess()->GetCurrentInteractionLength();
  if ( analogInteractionLength > DBL_MAX/10. ) return 0;

  // -- Analog cross-section is well-defined:
  G4double analogXS = 1./analogInteractionLength;

  // -- Choose a constant cross-section bias. But at this level, this factor can be made
  // -- direction dependent, like in the exponential transform MCNP case.
  G4double XStransformation = 2.0 ;
  
  G4BOptnChangeCrossSection* operation = fChangeCrossSectionOperations[callingProcess];
  // if ( ( callingProcess->GetPreviousOccurenceBiasingOperation() == operation ) &&
  //      ! operation->GetInteractionOccured() )
  //   {
  //     // -- operation proposed in previous step, but interaction did not occur
  //     // -- Update #int length for this step, and set new XS value:
  //     operation->UpdateForStep( callingProcess->GetPreviousStepSize() );
  //     operation->SetBiasedCrossSection( XStransformation * analogXS );
  //   }
  // else
  //   {
  //     // -- operation is new or interaction occured : we (re)sample:
  //     operation->SetBiasedCrossSection( XStransformation * analogXS );
  //     operation->Sample();
  //   }
  
  operation->SetBiasedCrossSection( XStransformation * analogXS );
  operation->Sample();

  // G4cout << callingProcess->GetProcessName() << " returning operation " << 
  // operation->GetName() << G4endl;

  return operation;
  
}


void GB01BOptrChangeCrossSection::
OperationApplied(const G4BiasingProcessInterface*           callingProcess, 
                 G4BiasingAppliedCase,
                 G4VBiasingOperation*             occurenceOperationApplied,
                 G4double,
                 G4VBiasingOperation*,    
                 const G4VParticleChange*                                  )
{
  G4BOptnChangeCrossSection* operation = fChangeCrossSectionOperations[callingProcess];
  if ( operation ==  occurenceOperationApplied ) operation->SetInteractionOccured();
}

