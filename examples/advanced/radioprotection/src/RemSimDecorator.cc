#include "RemSimVGeometryComponent.hh"
#include "RemSimDecorator.hh"

RemSimDecorator::RemSimDecorator(RemSimVGeometryComponent* geoComponent)
  : component(geoComponent)
{
}

RemSimDecorator::~RemSimDecorator()
{;}

void RemSimDecorator::ConstructComponent(G4VPhysicalVolume* motherVolume)
{
  //  component -> ConstructComponent(motherVolume) ; 
 }

void RemSimDecorator::DestroyComponent()
{
  //component -> DestroyComponent();
}
G4double  RemSimDecorator::GetDensity()
{;}

void  RemSimDecorator::ChangeMaterial(G4String)
{;}
