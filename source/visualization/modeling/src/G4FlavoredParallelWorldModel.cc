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
// $Id: G4FlavoredParallelWorldModel.cc,v 1.6 2002-11-20 17:21:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
    GetLogicalVolume () -> SetVisAttributes (G4VisAttributes::GetInvisible());
}

G4FlavoredParallelWorldModel::~G4FlavoredParallelWorldModel () {}
