//
// File name:     RadmonDetectorLayerVolumeItemUnion.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumeItemUnion.cc,v 1.1 2005-09-27 13:53:30 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorLayerVolumeItemUnion.hh"
#include "G4UnionSolid.hh"

G4VSolid *                                      RadmonDetectorLayerVolumeItemUnion :: Operate(G4VSolid * left, G4VSolid * right, G4RotationMatrix * relativeRotation, const G4ThreeVector & relativePosition)
{
 return new G4UnionSolid(left->GetName()+" + "+right->GetName(), left, right, relativeRotation, relativePosition);
}
