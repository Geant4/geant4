#ifndef RemSimDecorator_h
#define RemSimDecorator_h 1

#include "RemSimVGeometryComponent.hh"

class G4VPhysicalVolume;
class RemSimVGeometryComponent;

class RemSimDecorator: public RemSimVGeometryComponent
{
public:
  RemSimDecorator(RemSimVGeometryComponent*);
  ~RemSimDecorator();

  virtual void ConstructComponent(G4VPhysicalVolume*);
  virtual void DestroyComponent(); 
  virtual void ChangeThickness(G4double)=0;
  virtual void ChangeMaterial(G4String)=0;
  virtual void ChangeMother(G4VPhysicalVolume*)=0;
private:
   RemSimVGeometryComponent* component;
};
#endif
