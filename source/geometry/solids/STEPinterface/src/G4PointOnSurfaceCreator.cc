// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PointOnSurfaceCreator.cc,v 1.3 2000-02-25 16:36:19 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4PointOnSurfaceCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4PointOnSurfaceCreator.hh"
#include "G4GeometryTable.hh"

G4PointOnSurfaceCreator G4PointOnSurfaceCreator::csc;

G4PointOnSurfaceCreator::G4PointOnSurfaceCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4PointOnSurfaceCreator::~G4PointOnSurfaceCreator() {}

void G4PointOnSurfaceCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4Surface* srf=0;
  G4double pval1 = 0;
  G4double pval2 = 0;
  
  Ent.ResetAttributes();
  STEPattribute* Attr = Ent.NextAttribute();
  while(Attr->NonRefType() == STRING_TYPE ||
	Attr->NonRefType() == sdaiSTRING )
    Attr = Ent.NextAttribute();	

  // Get basis surface
  STEPentity* TmpEnt= *Attr->ptr.c;
  void * tmp =G4GeometryTable::CreateObject(*TmpEnt);
  srf = (G4Surface*)tmp;
  if (!tmp)
    G4cerr << "WARNING - G4PointOnSurfaceCreator::CreateG4Geometry" << G4endl
           << "\tUnexpected NULL pointer to G4Surface !" << G4endl;

  // Get parameter values
  Attr = Ent.NextAttribute();	
  pval1 = *Attr->ptr.r;
  Attr = Ent.NextAttribute();	
  pval2 = *Attr->ptr.r;

  G4double* dtmp = new G4double[2];
  dtmp[0] = pval1;
  dtmp[1] = pval2;  
  
  createdObject = dtmp;  
}

void G4PointOnSurfaceCreator::CreateSTEPGeometry(void* G4obj) {}
