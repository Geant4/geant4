// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FlavoredParallelWorldModel.cc,v 1.2 1999-01-11 00:48:51 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// P. Mora de Freitas et M.Verderi - 19 June 1998.
// Model for flavored parallel world volumes.

#include "G4FlavoredParallelWorldModel.hh"

#include "G4FlavoredParallelWorld.hh"
#include "G4VisAttributes.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

G4FlavoredParallelWorldModel::G4FlavoredParallelWorldModel
(G4FlavoredParallelWorld* FPW,
 G4int soughtDepth,
 const G4Transform3D& modelTransformation,
 const G4ModelingParameters* mp)
  : G4PhysicalVolumeModel (FPW -> GetThePhysicalVolumeWorld (),
			   soughtDepth,
			   modelTransformation,
			   mp),
    theFlavoredParallelWorld (FPW) 
{
  FPW -> GetThePhysicalVolumeWorld () ->
    GetLogicalVolume () -> SetVisAttributes (G4VisAttributes::Invisible);
}

G4FlavoredParallelWorldModel::~G4FlavoredParallelWorldModel () {}
