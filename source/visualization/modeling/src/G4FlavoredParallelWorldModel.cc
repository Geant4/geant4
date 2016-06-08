// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FlavoredParallelWorldModel.cc,v 1.3.8.1 1999/12/07 20:54:10 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// P. Mora de Freitas et M.Verderi - 19 June 1998.
// Model for flavored parallel world volumes.

#include "G4FlavoredParallelWorldModel.hh"

#include "G4VFlavoredParallelWorld.hh"
#include "G4VisAttributes.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

G4FlavoredParallelWorldModel::G4FlavoredParallelWorldModel
(G4VFlavoredParallelWorld* FPW,
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
