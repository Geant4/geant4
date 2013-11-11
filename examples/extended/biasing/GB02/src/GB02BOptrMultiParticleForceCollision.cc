#include "GB02BOptrMultiParticleForceCollision.hh"
#include "G4BiasingProcessInterface.hh"

#include "G4BOptrForceCollision.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

GB02BOptrMultiParticleForceCollision::GB02BOptrMultiParticleForceCollision()
  : G4VBiasingOperator("TestManyForceCollision")
{
}

void GB02BOptrMultiParticleForceCollision::AddParticle(G4String particleName)
{
  const G4ParticleDefinition* particle = 
    G4ParticleTable::GetParticleTable()->FindParticle( particleName );
  if ( particle == 0 ) 
    {
      G4cout << " ************* particle not found !!! : `" << particleName << "'" << G4endl;
      return;
    }
  G4BOptrForceCollision* optr = new G4BOptrForceCollision(particleName);
  fParticlesToBias.push_back( particle );
  fBOptrForParticle[ particle ] = optr;
}

G4VBiasingOperation* 
GB02BOptrMultiParticleForceCollision::
ProposeOccurenceBiasingOperation(const G4Track* track, 
                                 const G4BiasingProcessInterface* callingProcess)
{
  if ( fCurrentOperator ) return fCurrentOperator->
                            GetProposedOccurenceBiasingOperation(track, callingProcess);
  else                    return 0;
}

G4VBiasingOperation*
GB02BOptrMultiParticleForceCollision::
ProposeNonPhysicsBiasingOperation(const G4Track* track,
                                  const G4BiasingProcessInterface* callingProcess)
{
  if ( fCurrentOperator ) return fCurrentOperator->
                            GetProposedNonPhysicsBiasingOperation(track, callingProcess);
  else                    return 0;
}


void GB02BOptrMultiParticleForceCollision::StartTracking( const G4Track* track )
{
  const G4ParticleDefinition* definition = track->GetParticleDefinition();
  std::map < const G4ParticleDefinition*, G4BOptrForceCollision* > :: iterator
    it = fBOptrForParticle.find( definition );
  fCurrentOperator = 0;
  if ( it != fBOptrForParticle.end() ) fCurrentOperator = (*it).second;
}

void 
GB02BOptrMultiParticleForceCollision::
OperationApplied( const G4BiasingProcessInterface*         callingProcess,
                  G4BiasingAppliedCase                        biasingCase,
                  G4VBiasingOperation*                   operationApplied,
                  const G4VParticleChange*         particleChangeProduced )
{
  if ( fCurrentOperator ) fCurrentOperator->ReportOperationApplied( callingProcess,
                                                                    biasingCase,
                                                                    operationApplied,
                                                                    particleChangeProduced );
  
}

void
GB02BOptrMultiParticleForceCollision::
ExitBiasing( const G4Track*                           track,
             const G4BiasingProcessInterface* callingProcess )
{
  if ( fCurrentOperator ) fCurrentOperator->ExitingBiasing( track, callingProcess );
}


