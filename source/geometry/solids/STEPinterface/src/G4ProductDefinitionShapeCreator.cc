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
// $Id: G4ProductDefinitionShapeCreator.cc,v 1.5 2002-11-21 16:49:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4ProductDefinitionShapeCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4ProductDefinitionShapeCreator.hh"
#include "G4GeometryTable.hh"

G4ProductDefinitionShapeCreator G4ProductDefinitionShapeCreator::csc;

G4ProductDefinitionShapeCreator::G4ProductDefinitionShapeCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4ProductDefinitionShapeCreator::~G4ProductDefinitionShapeCreator() {}

G4ProductDefinitionShapeCreator G4ProductDefinitionShapeCreator::GetInstance()
{
  return csc;
}

void G4ProductDefinitionShapeCreator::CreateG4Geometry(STEPentity& Ent)
{
  // Made by L. Broglia
  STEPattribute *Attr;

  G4String attrName("description");
  Attr = GetNamedAttribute(attrName, Ent);

  // Get description
  SdaiString Tmpstring = *Attr->ptr.S; // ptr.S --> STRING_TYPE

  attrName = "definition";
  Attr = GetNamedAttribute(attrName, Ent);

  // Get definition
  SCLundefined Tmpdef = *Attr->ptr.u; // ptr.u --> UNKNOWN_TYPE

  // void *tmp = G4GeometryTable::CreateObject(*TmpEnt);
  // place = (G4Axis2Placement3D*)tmp;
}

void G4ProductDefinitionShapeCreator::CreateSTEPGeometry(void* G4obj)
{
}
