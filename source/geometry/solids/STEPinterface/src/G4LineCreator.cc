// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LineCreator.cc,v 1.2 2000-01-21 13:46:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4LineCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4LineCreator.hh"
#include "G4GeometryTable.hh"
#include "G4Line.hh"

G4LineCreator G4LineCreator::csc;

G4LineCreator::G4LineCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4LineCreator::~G4LineCreator() {}

void G4LineCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4ThreeVector origin;
  G4ThreeVector dir;
  
  Ent.ResetAttributes();
  STEPattribute *Attr;
  Attr = Ent.NextAttribute();    
  while (Attr->NonRefType() == STRING_TYPE)
    Attr = Ent.NextAttribute();

  // Get the line origin  
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  origin = *(G4ThreeVector*)tmp;

  // Get the line direction  
  Attr = Ent.NextAttribute();  
  TmpEnt = *Attr->ptr.c;
  tmp =G4GeometryTable::CreateObject(*TmpEnt);
  dir = *(G4ThreeVector*)tmp;


  G4Line* line = new G4Line();
  line->Init(origin, dir);
  createdObject = line;
}

void G4LineCreator::CreateSTEPGeometry(void* G4obj)
{
}
