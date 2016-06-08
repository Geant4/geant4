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
// Persistent class describing solid placements for boolean operations
//
// History:
// 10.11.99 Y.Morita, Initial creation

#include "G4PDisplacedSolid.hh"
#include "G4DisplacedSolid.hh"
#include "G4AffineTransform.hh"

#include "G4VSolid.hh"

G4PDisplacedSolid::G4PDisplacedSolid
                       ( HepRef(G4PVSolid) persCostituentSolid,
                         HepRef(G4PAffineTransform) pDirectTransform )
{
  fPtrSolid = persCostituentSolid;
  fDirectTransform = pDirectTransform;
}

G4PDisplacedSolid::~G4PDisplacedSolid()
{;}

HepRef(G4PVSolid) G4PDisplacedSolid::GetConstituentMovedSolid()
{ return fPtrSolid; }

G4VSolid* G4PDisplacedSolid::MakeTransientObject() const
{ return 0; }

G4VSolid* G4PDisplacedSolid::MakeTransientDisplacedSolid
                                     (G4VSolid* movedSolid) const
{
  G4AffineTransform aTransform = fDirectTransform->MakeTransientObject();

  G4VSolid* theSolid = new G4DisplacedSolid
                              ( GetName(), movedSolid, aTransform );
  return theSolid;
}


