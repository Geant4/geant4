// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OrientedEdgeCreator.cc,v 1.3 2000-02-25 16:36:19 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4OrientedEdgeCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4OrientedEdgeCreator.hh"
#include "G4GeometryTable.hh"

G4OrientedEdgeCreator G4OrientedEdgeCreator::csc;

G4OrientedEdgeCreator::G4OrientedEdgeCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4OrientedEdgeCreator::~G4OrientedEdgeCreator() {}

void G4OrientedEdgeCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4int orientation;
  G4Curve* crv=0;
  
  G4String attrName("edge_element");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);

  // Get curve
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  if (tmp)
  {
    crv = (G4Curve*)tmp; 

    // Get orientation
    Attr = Ent.NextAttribute();
    orientation = *Attr->ptr.i;// INTEGER_TYPE

    crv->SetSameSense(orientation);
  }
  else
    G4cerr << "WARNING - G4OrientedEdgeCreator::CreateG4Geometry" << G4endl
           << "\tUnexpected NULL pointer to G4Curve !" << G4endl
	   << "\tOriented Edge NOT created." << G4endl;

  createdObject = crv;
}

void G4OrientedEdgeCreator::CreateSTEPGeometry(void* G4obj)
{
}
