#ifndef RemSimVGeometryComponent_h
#define RemSimVGeometryComponent_h 1

#include "globals.hh"

class G4VPhysicalVolume; 
class RemSimVGeometryComponent
{
public:
  RemSimVGeometryComponent();
  virtual ~RemSimVGeometryComponent();

  virtual void ConstructComponent(G4VPhysicalVolume*)=0;
  virtual void DestroyComponent()=0;
};
#endif
