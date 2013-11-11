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
      G4cout << " ************* particle not found !!! : `" << particleName << "'" << G4endl;
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
  if ( fCurrentOperator ) return fCurrentOperator->
                            GetProposedOccurenceBiasingOperation(track, callingProcess);
  else                    return 0;
}


void GB01BOptrMultiParticleChangeCrossSection::StartTracking( const G4Track* track )
{
  const G4ParticleDefinition* definition = track->GetParticleDefinition();
  std::map < const G4ParticleDefinition*, GB01BOptrChangeCrossSection* > :: iterator
    it = fBOptrForParticle.find( definition );
  fCurrentOperator = 0;
  if ( it != fBOptrForParticle.end() ) fCurrentOperator = (*it).second;
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
  if ( fCurrentOperator ) fCurrentOperator->ReportOperationApplied( callingProcess,
                                                                    biasingCase,
                                                                    occurenceOperationApplied,
                                                                    weightForOccurenceInteraction,
                                                                    finalStateOperationApplied,
                                                                    particleChangeProduced );
}
