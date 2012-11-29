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
// $Id$
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
  fType = "G4FlavoredParallelWorldModel";
  FPW -> GetThePhysicalVolumeWorld () ->
    GetLogicalVolume () -> SetVisAttributes (G4VisAttributes::GetInvisible());
}

G4FlavoredParallelWorldModel::~G4FlavoredParallelWorldModel () {}
