// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParabolaCreator.cc,v 1.2 2000-01-21 13:46:04 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4ParabolaCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4ParabolaCreator.hh"
#include "G4GeometryTable.hh"

G4ParabolaCreator G4ParabolaCreator::csc;

G4ParabolaCreator::G4ParabolaCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4ParabolaCreator::~G4ParabolaCreator() {}

void G4ParabolaCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4double focal_dist;
  
  // L. Broglia
  // G4Placement* place;
  G4Axis2Placement3D* place;

  Ent.ResetAttributes();
  STEPattribute* Attr = Ent.NextAttribute();
  while(Attr->NonRefType() == STRING_TYPE ||
	Attr->NonRefType() == sdaiSTRING )
    Attr = Ent.NextAttribute();	

  // Get the placement
    STEPentity* TmpEnt= *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);

  // L. Broglia
  // place = (G4Placement*)tmp;
  place = (G4Axis2Placement3D*)tmp;
  
  // get semi axises
  Attr = Ent.NextAttribute();
  focal_dist = *Attr->ptr.r;

}

void G4ParabolaCreator::CreateSTEPGeometry(void* G4obj)
{
}
