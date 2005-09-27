//
// File name:     RadmonVDetectorLayerVolumeItemOperation.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorLayerVolumeItemOperation.cc,v 1.1 2005-09-27 13:53:30 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonVDetectorLayerVolumeItemOperation.hh"
#include "RadmonDetectorLayerVolumeItem.hh"
#include "G4VSolid.hh"

void                                            RadmonVDetectorLayerVolumeItemOperation :: Initialize(G4VSolid * solid, const G4RotationMatrix & rotation, const G4ThreeVector & position, G4bool ownSolid)
{
 Validate();

 leftRotation=rotation;
 invLeftRotation=rotation.inverse();
 leftPosition=position;
 leftSolid=solid;

 if (ownSolid)
  ownedSolids.push(solid);
}



void                                            RadmonVDetectorLayerVolumeItemOperation :: Initialize(G4VSolid * solid, const G4ThreeVector & position, G4bool ownSolid)
{
 Validate();

 leftRotation=G4RotationMatrix::IDENTITY;
 invLeftRotation=G4RotationMatrix::IDENTITY;
 leftPosition=position;
 leftSolid=solid;

 if (ownSolid)
  ownedSolids.push(solid);
}



void                                            RadmonVDetectorLayerVolumeItemOperation :: Initialize(G4VSolid * solid, G4bool ownSolid)
{
 Validate();

 leftRotation=G4RotationMatrix::IDENTITY;
 invLeftRotation=G4RotationMatrix::IDENTITY;
 leftPosition=G4ThreeVector(0, 0, 0);
 leftSolid=solid;

 if (ownSolid)
  ownedSolids.push(solid);
}



void                                            RadmonVDetectorLayerVolumeItemOperation :: Initialize(const RadmonDetectorLayerVolumeItem * item)
{
 Validate();

 leftSolid=item->GetSolid();

 Absolute(item, leftRotation, leftPosition);
 invLeftRotation=leftRotation.inverse();
}





G4VSolid *                                      RadmonVDetectorLayerVolumeItemOperation :: ApplyTo(G4VSolid * solid, const G4RotationMatrix & rotation, const G4ThreeVector & position, G4bool ownSolid, G4bool ownResult)
{
 if (ownSolid)
  ownedSolids.push(solid);

 G4ThreeVector relativePosition;
 ownedMatrices.push(G4RotationMatrix::IDENTITY);
 Merge(rotation, position, ownedMatrices.top(), relativePosition);

 G4VSolid * result(Operate(leftSolid, solid, &ownedMatrices.top(), relativePosition));
 
 if (ownResult)
  ownedSolids.push(result);

 return result;
}



void                                            RadmonVDetectorLayerVolumeItemOperation :: ApplyTo(RadmonDetectorLayerVolumeItem * item, G4bool ownResult)
{
 G4ThreeVector rightPosition;
 G4RotationMatrix rightRotation;
 Absolute(item, rightRotation, rightPosition);

 G4ThreeVector relativePosition;
 ownedMatrices.push(G4RotationMatrix::IDENTITY);
 Merge(rightRotation, rightPosition, ownedMatrices.top(), relativePosition);
 
 G4VSolid * result(Operate(leftSolid, item->GetSolid(), &ownedMatrices.top(), relativePosition));
 
 if (ownResult)
  ownedSolids.push(result);

 item->SetSolid(result);
 item->SetPosition(leftPosition);
 item->SetRotation(leftRotation);
}





                                                RadmonVDetectorLayerVolumeItemOperation :: RadmonVDetectorLayerVolumeItemOperation(G4VSolid * solid, const G4RotationMatrix & rotation, const G4ThreeVector & position, G4bool ownSolid)
:
 leftRotation(rotation),
 invLeftRotation(rotation.inverse()),
 leftPosition(position),
 leftSolid(solid)
{
 if (ownSolid)
  ownedSolids.push(solid);
}



                                                RadmonVDetectorLayerVolumeItemOperation :: RadmonVDetectorLayerVolumeItemOperation(G4VSolid * solid, const G4ThreeVector & position, G4bool ownSolid)
:
 leftRotation(G4RotationMatrix::IDENTITY),
 invLeftRotation(G4RotationMatrix::IDENTITY),
 leftPosition(position),
 leftSolid(solid)
{
 if (ownSolid)
  ownedSolids.push(solid);
}


                                                
                                                RadmonVDetectorLayerVolumeItemOperation :: RadmonVDetectorLayerVolumeItemOperation(G4VSolid * solid, G4bool ownSolid)
:
 leftRotation(G4RotationMatrix::IDENTITY),
 invLeftRotation(G4RotationMatrix::IDENTITY),
 leftPosition(0, 0, 0),
 leftSolid(solid)
{
 if (ownSolid)
  ownedSolids.push(solid);
}



                                                RadmonVDetectorLayerVolumeItemOperation :: RadmonVDetectorLayerVolumeItemOperation(const RadmonDetectorLayerVolumeItem * item)
:
 leftSolid(item->GetSolid())
{
 Absolute(item, leftRotation, leftPosition);
 invLeftRotation=leftRotation.inverse();
}



                                                RadmonVDetectorLayerVolumeItemOperation :: ~RadmonVDetectorLayerVolumeItemOperation()
{
 while (! ownedSolids.empty())
 {
  delete ownedSolids.top();
  ownedSolids.pop();
 }
}





inline void                                     RadmonVDetectorLayerVolumeItemOperation :: Validate()
{
 if (leftSolid!=0 || (!ownedSolids.empty()))
  G4Exception("RadmonVDetectorLayerVolumeItemOperation::Validate: Initialization happened twice.");
}





void                                            RadmonVDetectorLayerVolumeItemOperation :: Merge(const G4RotationMatrix & rightRotation, const G4ThreeVector & rightPosition, G4RotationMatrix & relativeRotation, G4ThreeVector & relativePosition) const
{
 // X_a = O_r + R_r * X_r                               X_r --> X_a
 // X_a = O_l + R_l * X_l                               X_l --> X_a
 // X_l = R_l^-1 * (O_r - O_l) + R_l^-1 * R_r * X_r     X_r --> X_l
 
 relativeRotation=rightRotation;
 relativeRotation.transform(invLeftRotation);
 
 relativePosition=rightPosition;
 relativePosition-=leftPosition;
 relativePosition.transform(invLeftRotation);
}



void                                            RadmonVDetectorLayerVolumeItemOperation :: Absolute(const RadmonDetectorLayerVolumeItem * item, G4RotationMatrix & rotation, G4ThreeVector & position) const
{
 rotation=item->GetRotation();
 position=item->GetPosition();

 for(;;)
 {
  item=item->GetMotherVolumeItem();
  if (!item)
   return;
   
  rotation.transform(item->GetRotation());
  position.transform(item->GetRotation());
 }
}

