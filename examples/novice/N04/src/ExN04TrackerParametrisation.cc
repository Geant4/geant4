
#include "ExN04TrackerParametrisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"

ExN04TrackerParametrisation::ExN04TrackerParametrisation()
{

#include "ExN04DetectorParameterDef.icc"

}

ExN04TrackerParametrisation::~ExN04TrackerParametrisation()
{;}

void ExN04TrackerParametrisation::ComputeTransformation
(const G4int copyNo,G4VPhysicalVolume *physVol) const
{
  G4ThreeVector origin;
  physVol->SetTranslation(origin);
}

void ExN04TrackerParametrisation::ComputeDimensions
(G4Tubs & trackerLayer, const G4int copyNo,
 const G4VPhysicalVolume * physVol) const
{
  trackerLayer.SetInnerRadius(tracker_radius[copyNo]);
  trackerLayer.SetOuterRadius(tracker_radius[copyNo]+tracker_thick);
  trackerLayer.SetZHalfLength(tracker_length[copyNo]);
  trackerLayer.SetStartPhiAngle(trkTubs_sphi);
  trackerLayer.SetDeltaPhiAngle(trkTubs_dphi);
}
