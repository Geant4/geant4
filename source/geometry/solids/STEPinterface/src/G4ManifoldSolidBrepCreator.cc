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
// $Id: G4ManifoldSolidBrepCreator.cc,v 1.5 2002-11-21 16:49:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4ManifoldSolidBrepCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4ManifoldSolidBrepCreator.hh"
#include "G4GeometryTable.hh"
#include "G4BREPSolid.hh"

G4ManifoldSolidBrepCreator G4ManifoldSolidBrepCreator::csc;

G4ManifoldSolidBrepCreator::G4ManifoldSolidBrepCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4ManifoldSolidBrepCreator::~G4ManifoldSolidBrepCreator() {}

G4ManifoldSolidBrepCreator G4ManifoldSolidBrepCreator::GetInstance()
{
  return csc;
}

void G4ManifoldSolidBrepCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4BREPSolid* sld;
  
  G4String attrName("outer");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);

  // Get solid
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  if (!tmp)
    G4cerr << "WARNING - G4ManifoldSolidBrepCreator::CreateG4Geometry" << G4endl
           << "\tUnexpected NULL pointer to BREP Solid !" << G4endl
	   << "\tBREP Solid NOT created." << G4endl;

  sld = (G4BREPSolid*)tmp; 
  createdObject = sld;
}

void G4ManifoldSolidBrepCreator::CreateSTEPGeometry(void* G4obj)
{
}
