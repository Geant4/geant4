#ifndef exGPSGeometryConstruction1_h
#define exGPSGeometryConstruction1_h 1
//
#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;

////////////////////////////////////////////////////////////////////////////////
//
class exGPSGeometryConstruction : public G4VUserDetectorConstruction
{
  public:
    exGPSGeometryConstruction ();
    ~exGPSGeometryConstruction ();

  public:
     G4VPhysicalVolume *Construct ();

  private:

  G4VPhysicalVolume * universe_phys ;
  G4VPhysicalVolume * Al_phys ;
  G4VPhysicalVolume * aSphere_phys;

};
////////////////////////////////////////////////////////////////////////////////
#endif

