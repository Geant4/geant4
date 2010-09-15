#include "G4InclUtils.hh"
#include "G4NucleiProperties.hh"

G4double G4InclUtils::calculate4MomentumScaling(G4int A, G4int Z, G4double excitationE, G4double kineticE,
						G4double px, G4double py, G4double pz)
{
  G4double nuclearMass = (G4NucleiProperties::GetNuclearMass(A, Z) / MeV + excitationE) * MeV;
  return std::sqrt(kineticE*kineticE + 2.0 * kineticE * nuclearMass)/std::sqrt(px*px + py*py + pz*pz);
}
