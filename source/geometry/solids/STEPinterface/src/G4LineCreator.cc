// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LineCreator.cc,v 1.3 2000-02-25 16:36:19 gcosmo Exp $
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
  if (tmp)
    origin = *(G4ThreeVector*)tmp;
  else
    G4cerr << "WARNING - G4LineCreator::CreateG4Geometry" << G4endl
           << "\tUnexpected NULL origin (G4ThreeVector) for G4Line !" << G4endl
	   << "\tG4Line origin set to (0,0,0)." << G4endl;

  // Get the line direction  
  Attr = Ent.NextAttribute();  
  TmpEnt = *Attr->ptr.c;
  tmp =G4GeometryTable::CreateObject(*TmpEnt);
  if (tmp)
    dir = *(G4ThreeVector*)tmp;
  else
    G4cerr << "WARNING - G4LineCreator::CreateG4Geometry" << G4endl
           << "\tUnexpected NULL direction (G4ThreeVector) for G4Line !" << G4endl
	   << "\tG4Line direction set to (0,0,0)." << G4endl;

  G4Line* line = new G4Line();
  line->Init(origin, dir);
  createdObject = line;
}

void G4LineCreator::CreateSTEPGeometry(void* G4obj)
{
}
