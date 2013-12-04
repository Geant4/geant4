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
  
  if ( fParticleToBias == 0 )
    {
      G4ExceptionDescription ed;
      ed << "Particle `" << particleName << "' not found !" << G4endl;
      G4Exception("GB01BOptrChangeCrossSection(...)",
                  "exGB01.01",
                  JustWarning,
                  ed);
    }
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

  // -----------------------------------------------------------------------
  // -- Setup stage ( Might consider moving this in a less called method. ):
  // -----------------------------------------------------------------------
  // -- Start by collecting processes under biasing, create needed biasing operations
  // -- and assocation operations to these processes:
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



  // ---------------------------------------------------------------------
  // -- select and setup the biasing operation for current callingProcess:
  // ---------------------------------------------------------------------
  // -- Check if the analog cross-section well defined : for example, the conversion
  // -- process for a gamma below e+e- creation threshold has an DBL_MAX interaction
  // -- length. Nothing is done in this case (ie, let analog process to deal with the case)
  G4double analogInteractionLength =  
    callingProcess->GetWrappedProcess()->GetCurrentInteractionLength();
  if ( analogInteractionLength > DBL_MAX/10. ) return 0;

  // -- Analog cross-section is well-defined:
  G4double analogXS = 1./analogInteractionLength;

  // -- Choose a constant cross-section bias. But at this level, this factor can be made
  // -- direction dependent, like in the exponential transform MCNP case, or it
  // -- can be chosen differently, depending on the process, etc.
  G4double XStransformation = 2.0 ;
  
  // -- fetch the operation associated to this callingProcess:
  G4BOptnChangeCrossSection*   operation = fChangeCrossSectionOperations[callingProcess];
  // -- get the operation that was proposed to the process in the previous step:
  G4VBiasingOperation* previousOperation = callingProcess->GetPreviousOccurenceBiasingOperation();
  
  // -- now setup the operation to be returned to the process: this
  // -- consists in setting the biased cross-section, and in asking
  // -- the operation to sample its exponential interaction law.
  // -- To do this, to first order, the two lines:
  //        operation->SetBiasedCrossSection( XStransformation * analogXS );
  //        operation->Sample();
  // -- are correct and sufficient.
  // -- But, to avoid having to shoot a random number each time, we sample
  // -- only on the first time the operation is proposed, or if the interaction
  // -- occured. If the interaction did not occur for the process in the previous,
  // -- we update the number of interaction length instead of resampling.
  if ( previousOperation == 0 )
    {
      operation->SetBiasedCrossSection( XStransformation * analogXS );
      operation->Sample();
    }
  else
    {
      if (  previousOperation != operation )
        {
          // -- should not happen !
          G4ExceptionDescription ed;
          ed << " Logic problem in operation handling !" << G4endl;
          G4Exception("GB01BOptrChangeCrossSection::ProposeOccurenceBiasingOperation(...)",
                      "exGB01.02",
                      JustWarning,
                      ed);
          return 0;
        }
      if ( operation->GetInteractionOccured() )
        {
          operation->SetBiasedCrossSection( XStransformation * analogXS );
          operation->Sample();
        }
      else
        {
          operation->UpdateForStep( callingProcess->GetPreviousStepSize() );
          operation->SetBiasedCrossSection( XStransformation * analogXS );
        }
    }
  
  //  operation->SetBiasedCrossSection( XStransformation * analogXS );
  //  operation->Sample();

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

