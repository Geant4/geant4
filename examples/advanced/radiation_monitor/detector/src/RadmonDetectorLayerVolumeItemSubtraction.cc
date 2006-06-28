//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// File name:     RadmonDetectorLayerVolumeItemSubtraction.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumeItemSubtraction.cc,v 1.2 2006-06-28 13:50:32 gunter Exp $
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
