#include "G4InclUtils.hh"
#include "G4NucleiProperties.hh"

G4double G4InclUtils::calculate4MomentumScaling(G4int A, G4int Z, G4double excitationE, G4double kineticE,
						G4double px, G4double py, G4double pz)
{
  G4double nuclearMass = (G4NucleiProperties::GetNuclearMass(A, Z) / MeV + excitationE) * MeV;
  G4double p2 = px*px + py*py + pz*pz;
  if(p2 > 0.0)
    return std::sqrt(kineticE*kineticE + 2.0 * kineticE * nuclearMass)/std::sqrt(p2);
  else // if p2 <= 0.0 we have incorrect input. Returning 1.0 allows us to attempt to recover from the error.
    return 1.0;
}
