// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ShapeDefinitionRepresentationCreator.cc,v 1.4 2000-02-25 16:36:20 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4ShapeDefinitionRepresentationCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4ShapeDefinitionRepresentationCreator.hh"
#include "G4GeometryTable.hh"

G4ShapeDefinitionRepresentationCreator
  G4ShapeDefinitionRepresentationCreator::csc;

G4ShapeDefinitionRepresentationCreator::
  G4ShapeDefinitionRepresentationCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4ShapeDefinitionRepresentationCreator::
  ~G4ShapeDefinitionRepresentationCreator() {}

void G4ShapeDefinitionRepresentationCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4String attrName("used_representation");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  G4PlacedSolid* ps = (G4PlacedSolid*)tmp;
  if (!tmp)
    G4cerr << "WARNING - G4ShapeDefinitionRepresentationCreator::CreateG4Geometry" << G4endl
           << "\tUnexpected NULL pointer to G4PlacedSolid !" << G4endl
	   << "\tEntity NOT created." << G4endl;

  createdObject = ps;

/*
  // Output messages
  G4cout<<"Shape definition representation created with "
	<<"the Advanced Brep shape :"<<G4endl;

  G4cout<<" - BBox of the solid :"<<G4endl;
  G4cout<< "    box min: "
	<< ((G4BREPSolid*)ps->GetSolid())->GetBBox()->GetBoxMin().x()<< " "
	<< ((G4BREPSolid*)ps->GetSolid())->GetBBox()->GetBoxMin().y()<< " "
	<< ((G4BREPSolid*)ps->GetSolid())->GetBBox()->GetBoxMin().z()<< G4endl;

  G4cout<< "    box max: "
	<< ((G4BREPSolid*)ps->GetSolid())->GetBBox()->GetBoxMax().x()<< " "
	<< ((G4BREPSolid*)ps->GetSolid())->GetBBox()->GetBoxMax().y()<< " "
	<< ((G4BREPSolid*)ps->GetSolid())->GetBBox()->GetBoxMax().z()<< G4endl;

*/
}

void G4ShapeDefinitionRepresentationCreator::CreateSTEPGeometry(void* G4obj)
{
}
