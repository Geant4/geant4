#ifndef RemSimDecorator_h
#define RemSimDecorator_h 1

#include "RemSimVGeometryComponent.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class RemSimVGeometryComponent;

class RemSimDecorator: public RemSimVGeometryComponent
{
public:
  RemSimDecorator(RemSimVGeometryComponent*);
    ~RemSimDecorator();

  virtual void ConstructComponent(G4VPhysicalVolume*);
  virtual void DestroyComponent(); 
  virtual G4double GetDensity();
  virtual void ChangeMaterial(G4String);
private:
   RemSimVGeometryComponent* component;
};
#endif
