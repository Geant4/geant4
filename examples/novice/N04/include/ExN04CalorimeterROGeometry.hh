#ifndef ExN04CalorimeterROGeometry_h
#define ExN04CalorimeterROGeometry_h 1

#include "G4VReadOutGeometry.hh"

class ExN04CalorimeterROGeometry : public G4VReadOutGeometry
{
public:
  ExN04CalorimeterROGeometry();
  ExN04CalorimeterROGeometry(G4String);
  ~ExN04CalorimeterROGeometry();

private:
  G4VPhysicalVolume* Build();

#include "ExN04DetectorParameterDef.hh"

};

#endif
