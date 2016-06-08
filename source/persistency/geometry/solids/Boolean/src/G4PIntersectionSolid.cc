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
// Persistent class for description of intersection of two CSG solids
//
// History: 
// 10.11.99 Y.Morita, Initial creation

#include "G4PIntersectionSolid.hh"
#include "G4IntersectionSolid.hh"


G4PIntersectionSolid::G4PIntersectionSolid
                            ( const G4String& pName,
                              HepRef(G4PVSolid) persSolidA,
                              HepRef(G4PVSolid) persSolidB )
 : G4PBooleanSolid(pName, persSolidA, persSolidB)
{;}

G4PIntersectionSolid::~G4PIntersectionSolid()
{;}

G4VSolid* G4PIntersectionSolid::MakeTransientObject() const
{ return 0; }

G4VSolid* G4PIntersectionSolid::MakeTransientBooleanSolid(
                              G4VSolid* aSolidA,
                              G4VSolid* aSolidB ) const
{
  return new G4IntersectionSolid( GetName(), aSolidA, aSolidB );
}


