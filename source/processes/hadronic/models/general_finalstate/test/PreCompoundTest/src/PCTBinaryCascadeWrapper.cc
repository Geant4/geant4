#include "PCTBinaryCascadeWrapper.hh"

#include "G4Track.hh"
#include "G4Nucleus.hh"
#include "G4DynamicParticle.hh"
#include "G4Step.hh"
#include "G4HadFinalState.hh"

G4ReactionProductVector * PCTBinaryCascadeWrapper::DeExcite(const PCTProjectile * theProjectile,
							    const G4int targetA, const G4int targetZ)
{
  // Build the Target Nucleus
  G4Nucleus aNucleus(targetA, targetZ);

  // Build the Projectile Track
  G4DynamicParticle projDynamicParticle(theProjectile->GetDefinition(), theProjectile->GetMomentum());
  G4double time(0.0);
  G4ThreeVector position(0.0,0.0,0.0);
  G4Track projTrack(&projDynamicParticle, time, position);
  G4Step * aStep = new G4Step();
  projTrack.SetStep(aStep);

  // Produce secondaries
  G4HadFinalState * thePChResult = theCascade.ApplyYourself(projTrack, aNucleus);

  // Trasform Back G4VParticleChange to G4ReactionProduct
  G4ReactionProductVector * result = new G4ReactionProductVector();
  result->reserve(thePChResult->GetNumberOfSecondaries());

  for (G4int i = 0; i < thePChResult->GetNumberOfSecondaries(); i++)
    {
      G4HadSecondary * secondary = thePChResult->GetSecondary(i);
      G4ReactionProduct * arp = new G4ReactionProduct(secondary->GetParticle()->GetDefinition());
      arp->SetMass(secondary->GetParticle()->GetMass());
      arp->SetMomentum(secondary->GetParticle()->GetMomentum());
      arp->SetKineticEnergy(secondary->GetParticle()->GetKineticEnergy());
      result->push_back(arp);
    }
  delete thePChResult;
  delete aStep;
  return result;
}
