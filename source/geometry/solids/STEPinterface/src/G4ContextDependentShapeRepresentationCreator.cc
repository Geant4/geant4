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
//
// $Id: G4ContextDependentShapeRepresentationCreator.cc,v 1.4 2001-07-11 10:00:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ContextDependentShapeRepresentationCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4ContextDependentShapeRepresentationCreator.hh"
#include "G4GeometryTable.hh"

G4ContextDependentShapeRepresentationCreator
  G4ContextDependentShapeRepresentationCreator::csc;

G4ContextDependentShapeRepresentationCreator::
  G4ContextDependentShapeRepresentationCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4ContextDependentShapeRepresentationCreator::
  ~G4ContextDependentShapeRepresentationCreator() {}

void G4ContextDependentShapeRepresentationCreator::
  CreateG4Geometry(STEPentity& Ent)
{
  G4String attrName("representation_relation");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  if (!tmp)
    G4cerr << "WARNING - G4ContextDependentShapeRepresentationCreator::CreateG4Geometry" << G4endl
           << "\tNULL G4PlacedSolidVector returned." << G4endl;

  G4PlacedSolidVector* psV = (G4PlacedSolidVector*)tmp;
    
  // L. Broglia
  // For the moment, the represented_product_relation is not needed
  // So don`t create it
  attrName="represented_product_relation";
  Attr = GetNamedAttribute(attrName, Ent);
  TmpEnt = *Attr->ptr.c;
  // tmp =G4GeometryTable::CreateObject(*TmpEnt);

  createdObject = psV;
}

void G4ContextDependentShapeRepresentationCreator::
  CreateSTEPGeometry(void* obj)
{
}
