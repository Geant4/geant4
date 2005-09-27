//
// File name:     RadmonDetectorLayerVolumeItemSubtraction.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumeItemSubtraction.cc,v 1.1 2005-09-27 13:53:30 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorLayerVolumeItemSubtraction.hh"
#include "G4SubtractionSolid.hh"
#include "G4DisplacedSolid.hh"

G4VSolid *                                      RadmonDetectorLayerVolumeItemSubtraction :: Operate(G4VSolid * left, G4VSolid * right, G4RotationMatrix * relativeRotation, const G4ThreeVector & relativePosition)
{
 if (opMode==leftMinusRight)
  return new G4SubtractionSolid(left->GetName()+" - "+right->GetName(), left, right, relativeRotation, relativePosition);

 G4RotationMatrix * inverseRelativeRotation(AllocateMatrix());
 (*inverseRelativeRotation)=(*relativeRotation);
 inverseRelativeRotation->invert();
 G4VSolid *solid(new G4SubtractionSolid(right->GetName()+" - "+left->GetName(), right, left, inverseRelativeRotation, -relativePosition));
 OwnSolid(solid);
 
 return new G4DisplacedSolid("- "+left->GetName()+" + "+right->GetName(), solid, relativeRotation, relativePosition);
}
