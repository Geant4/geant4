// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GeometryCreator.cc,v 1.5 2000-11-10 17:44:46 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4GeometryCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include <instmgr.h>
#include <STEPcomplex.h>

#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

InstMgr G4GeometryCreator::instanceManager;
G4int G4GeometryCreator::objectId = 0;

G4GeometryCreator::G4GeometryCreator() {}
G4GeometryCreator::~G4GeometryCreator() {}

G4bool G4GeometryCreator::operator==(const G4GeometryCreator&)
{
  return 0;
}

void* G4GeometryCreator::GetCreatedObject()
{
  return createdObject;
}

InstMgr* G4GeometryCreator::GetInstanceManager() const
{
  return &instanceManager;
}

STEPattribute*
G4GeometryCreator::GetNamedAttribute(const G4String& attrName, STEPentity& Ent)
{
  STEPattribute* attr;

  Ent.ResetAttributes();
  G4int attrCount = Ent.AttributeCount();

  for(G4int a=0;a<attrCount;a++)
    {
      attr = Ent.NextAttribute();
      if(attr->Name() == attrName)
	return attr;
    }
      
  G4String err = "\nCannot find attribute " + G4String(attrName)
                 + " in entity " + Ent.EntityName();
  G4Exception(err);
  return 0; // NULL
}

STEPentity*
G4GeometryCreator::GetNamedEntity(const G4String& entName, STEPentity& Ent)
{
  // Ent is a complex entity

  if(Ent.IsComplex())
    {
      Ent.ResetAttributes();
      STEPcomplex *complexEnt = (STEPcomplex*)&Ent;
      STEPentity* subEnt=0;
      
      G4int entCount = complexEnt->AttributeCount();

      for(G4int a=0;a<entCount;a++)
	{
	  subEnt = complexEnt->sc;
	  if(subEnt->EntityName() == entName)
	    return subEnt;
	}
     } 
  G4String err = "\nCannot find entity " + G4String(entName)
                 + " in complex entity " + Ent.EntityName();
  G4Exception(err);

  return 0;    
}

/*
G4int[][] G4GeometryCreator::GetIds(GenericAggregate& aggr)
{
  // Gets the file IDs of entities from the string.
  // hack for 2D-aggregates which do not work properly in the NIST toolkit
    GenericAggrNode* gNode = (GenericAggrNode*)aggr->GetHead();

  SCLstring str;
  const char* aggrStr = gNode->asStr(str); 
}
*/
