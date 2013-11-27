#include "GB01BOptrMultiParticleChangeCrossSection.hh"
#include "G4BiasingProcessInterface.hh"

#include "GB01BOptrChangeCrossSection.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

GB01BOptrMultiParticleChangeCrossSection::GB01BOptrMultiParticleChangeCrossSection()
  : G4VBiasingOperator("TestManyExponentialTransform")
{}

void GB01BOptrMultiParticleChangeCrossSection::AddParticle(G4String particleName)
{
  const G4ParticleDefinition* particle =
    G4ParticleTable::GetParticleTable()->FindParticle( particleName );
  
  if ( particle == 0 )
    {
      G4ExceptionDescription ed;
      ed << "Particle `" << particleName << "' not found !" << G4endl;
      G4Exception("GB01BOptrMultiParticleChangeCrossSection::AddParticle(...)",
                  "exGB01.02",
                  JustWarning,
                  ed);
      return;
    }
  
  GB01BOptrChangeCrossSection* optr = new GB01BOptrChangeCrossSection(particleName);
  fParticlesToBias.push_back( particle );
  fBOptrForParticle[ particle ] = optr;
}

G4VBiasingOperation*
GB01BOptrMultiParticleChangeCrossSection::
ProposeOccurenceBiasingOperation(const G4Track* track,
                                 const G4BiasingProcessInterface* callingProcess)
{
  // -- examples of limitations imposed to applied the biasing:
  // -- limit application of biasing to primary particles only:
  if ( track->GetParentID() != 0 ) return 0;
  // -- limit to at most 5 biased interactions:
  if ( fnInteractions > 4 )        return 0;
  // -- and limit to a weight of at least 0.05:
  if ( track->GetWeight() < 0.05 ) return 0;

  if ( fCurrentOperator ) return fCurrentOperator->
                            GetProposedOccurenceBiasingOperation(track, callingProcess);
  else                    return 0;
}


void GB01BOptrMultiParticleChangeCrossSection::StartTracking( const G4Track* track )
{
  // -- fetch the underneath biasing operator, if any, for the current particle type:
  const G4ParticleDefinition* definition = track->GetParticleDefinition();
  std::map < const G4ParticleDefinition*, GB01BOptrChangeCrossSection* > :: iterator
    it = fBOptrForParticle.find( definition );
  fCurrentOperator = 0;
  if ( it != fBOptrForParticle.end() ) fCurrentOperator = (*it).second;

  // -- reset count for number of biased interactions:
  fnInteractions = 0;
}

void 
GB01BOptrMultiParticleChangeCrossSection::
OperationApplied( const G4BiasingProcessInterface*               callingProcess, 
                  G4BiasingAppliedCase                              biasingCase,
                  G4VBiasingOperation*                occurenceOperationApplied, 
                  G4double                        weightForOccurenceInteraction,
                  G4VBiasingOperation*               finalStateOperationApplied, 
                  const G4VParticleChange*               particleChangeProduced )
{
  // -- count number of biased interactions:
  fnInteractions++;

  // -- inform the underneath biasing operator that a biased interaction occured:
  if ( fCurrentOperator ) fCurrentOperator->ReportOperationApplied( callingProcess,
                                                                    biasingCase,
                                                                    occurenceOperationApplied,
                                                                    weightForOccurenceInteraction,
                                                                    finalStateOperationApplied,
                                                                    particleChangeProduced );
}
