#include "G4GDMLParameterisation.hh"

G4int G4GDMLParameterisation::getSize() const {

   return (G4int)parameterList.size();
}

void G4GDMLParameterisation::addParameter(const PARAMETER& newParameter) {

   parameterList.push_back(newParameter);
}

void G4GDMLParameterisation::ComputeTransformation(const G4int index,G4VPhysicalVolume* physvol) const {

   physvol->SetTranslation(parameterList[index].position);
   physvol->SetRotation(parameterList[index].pRot);
}

void G4GDMLParameterisation::ComputeDimensions(G4Box& box,const G4int index,const G4VPhysicalVolume*) const {

  box.SetXHalfLength(parameterList[index].dimension[0]);
  box.SetYHalfLength(parameterList[index].dimension[1]);
  box.SetZHalfLength(parameterList[index].dimension[2]);
}
