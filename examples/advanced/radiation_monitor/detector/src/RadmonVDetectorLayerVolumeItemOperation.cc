//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// File name:     RadmonVDetectorLayerVolumeItemOperation.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorLayerVolumeItemOperation.cc,v 1.3 2006/06/29 16:14:15 gunter Exp $
// Tag:           $Name: geant4-08-01 $
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

