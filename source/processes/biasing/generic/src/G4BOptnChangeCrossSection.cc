#include "G4BOptnChangeCrossSection.hh"
#include "G4InteractionLawPhysical.hh"



G4BOptnChangeCrossSection::G4BOptnChangeCrossSection(G4String name)
  : G4VBiasingOperation(name)
{
  fBiasedExponentialLaw = new G4InteractionLawPhysical("LawForOperation"+name);
}

G4BOptnChangeCrossSection::~G4BOptnChangeCrossSection()
{}

const G4VBiasingInteractionLaw* G4BOptnChangeCrossSection::ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface* )
{
  return fBiasedExponentialLaw;
}

void G4BOptnChangeCrossSection::SetBiasedCrossSection(G4double xst)
{
  fBiasedExponentialLaw->SetPhysicalCrossSection( xst );
}

G4double G4BOptnChangeCrossSection::GetBiasedCrossSection() const
{
  return fBiasedExponentialLaw->GetPhysicalCrossSection();
}

void G4BOptnChangeCrossSection::Sample()
{
  fInteractionOccured = false;
  fBiasedExponentialLaw->Sample();
}

void G4BOptnChangeCrossSection::UpdateForStep( G4double truePathLength )
{
  fBiasedExponentialLaw->UpdateForStep( truePathLength );
}
