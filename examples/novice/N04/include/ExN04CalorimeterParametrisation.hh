
#ifndef ExN04CalorimeterParametrisation_H
#define ExN04CalorimeterParametrisation_H 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
class G4VPhysicalVolume;
class G4Tubs;

class ExN04CalorimeterParametrisation : public G4VPVParameterisation
{ 
  public:
    ExN04CalorimeterParametrisation();
    ~ExN04CalorimeterParametrisation();
    void ComputeTransformation
    (const G4int copyNo,G4VPhysicalVolume *physVol) const;
    void ComputeDimensions
    (G4Tubs & calorimeterLayer, const G4int copyNo,
      const G4VPhysicalVolume * physVol) const;

  private:

#include "ExN04DetectorParameterDef.hh"

};

#endif


