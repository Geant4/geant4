#include "G4BOptnCloning.hh"


G4BOptnCloning::G4BOptnCloning(G4String name)
  : G4VBiasingOperation(name),
    fParticleChange()
{}

G4BOptnCloning::~G4BOptnCloning()
{}

G4VParticleChange*  G4BOptnCloning::GenerateBiasingFinalState( const G4Track* track,
								    const G4Step*       )
{
  fParticleChange.Initialize(*track);
  fParticleChange.ProposeParentWeight( fClone1W );
  fParticleChange.SetSecondaryWeightByProcess(true);
  fParticleChange.SetNumberOfSecondaries(1);
  G4Track* clone = new G4Track( *track );
  clone->SetWeight( fClone2W );
  fParticleChange.AddSecondary( clone );
  return &fParticleChange;
}
