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

void G4GDMLParameterisation::ComputeDimensions(G4Trd& trd,const G4int index,const G4VPhysicalVolume*) const {

   trd.SetXHalfLength1(parameterList[index].dimension[0]);
   trd.SetXHalfLength2(parameterList[index].dimension[1]);
   trd.SetYHalfLength1(parameterList[index].dimension[2]);
   trd.SetYHalfLength2(parameterList[index].dimension[3]);
   trd.SetZHalfLength(parameterList[index].dimension[4]);
}

void G4GDMLParameterisation::ComputeDimensions(G4Trap& trap,const G4int index,const G4VPhysicalVolume*) const {

   trap.SetAllParameters(parameterList[index].dimension[0], // Dz
                         parameterList[index].dimension[1], // Theta
                         parameterList[index].dimension[2], // Phi
                         parameterList[index].dimension[3], // Dy1
                         parameterList[index].dimension[4], // Dx1
                         parameterList[index].dimension[5], // Dx2
                         parameterList[index].dimension[6], // pAlp1,
                         parameterList[index].dimension[7], // pDy2,
                         parameterList[index].dimension[8], // pDx3,
                         parameterList[index].dimension[9], // pDx4,
                         parameterList[index].dimension[10]); // pAlp2
}

void G4GDMLParameterisation::ComputeDimensions(G4Tubs& tubs,const G4int index,const G4VPhysicalVolume*) const {

   tubs.SetInnerRadius(parameterList[index].dimension[0]);
   tubs.SetOuterRadius(parameterList[index].dimension[1]);
   tubs.SetZHalfLength(parameterList[index].dimension[2]);
   tubs.SetStartPhiAngle(parameterList[index].dimension[3]);
   tubs.SetDeltaPhiAngle(parameterList[index].dimension[4]);
}

void G4GDMLParameterisation::ComputeDimensions(G4Cons& cons,const G4int index,const G4VPhysicalVolume*) const {

   cons.SetInnerRadiusMinusZ(parameterList[index].dimension[0]);
   cons.SetOuterRadiusMinusZ(parameterList[index].dimension[1]);
   cons.SetInnerRadiusPlusZ(parameterList[index].dimension[2]);
   cons.SetOuterRadiusPlusZ(parameterList[index].dimension[3]);
   cons.SetZHalfLength(parameterList[index].dimension[4]);
   cons.SetStartPhiAngle(parameterList[index].dimension[5]);
   cons.SetDeltaPhiAngle(parameterList[index].dimension[6]);
}


