
#ifndef ExN04TrackerParametrisation_H
#define ExN04TrackerParametrisation_H 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
class G4VPhysicalVolume;
class G4Tubs;

class ExN04TrackerParametrisation : public G4VPVParameterisation
{ 
  public:
    ExN04TrackerParametrisation();
    ~ExN04TrackerParametrisation();
    void ComputeTransformation
    (const G4int copyNo,G4VPhysicalVolume *physVol) const;
    void ComputeDimensions
    (G4Tubs & trackerLayer, const G4int copyNo,
      const G4VPhysicalVolume * physVol) const;

  private:

#include "ExN04DetectorParameterDef.hh"

};

#endif


