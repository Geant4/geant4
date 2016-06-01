
#include "ExN04CalorimeterParametrisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"

ExN04CalorimeterParametrisation::ExN04CalorimeterParametrisation()
{

#include "ExN04DetectorParameterDef.icc"

}

ExN04CalorimeterParametrisation::~ExN04CalorimeterParametrisation()
{;}

void ExN04CalorimeterParametrisation::ComputeTransformation
(const G4int copyNo,G4VPhysicalVolume *physVol) const
{
  G4ThreeVector origin;
  physVol->SetTranslation(origin);
}

void ExN04CalorimeterParametrisation::ComputeDimensions
(G4Tubs & calorimeterLayer, const G4int copyNo,
 const G4VPhysicalVolume * physVol) const
{
  G4double innerRad = caloTubs_rmin
              + copyNo*(absorber_thick+scinti_thick);
  calorimeterLayer.SetInnerRadius(innerRad);
  calorimeterLayer.SetOuterRadius(innerRad+absorber_thick);
  calorimeterLayer.SetZHalfLength(caloTubs_dz);
  calorimeterLayer.SetStartPhiAngle(caloTubs_sphi);
  calorimeterLayer.SetDeltaPhiAngle(caloTubs_dphi);
}
