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
// $Id: G4EllipseCreator.cc,v 1.5 2002-11-21 16:49:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4EllipseCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4EllipseCreator.hh"
#include "G4GeometryTable.hh"

G4EllipseCreator G4EllipseCreator::csc;

G4EllipseCreator::G4EllipseCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4EllipseCreator::~G4EllipseCreator() {}

G4EllipseCreator G4EllipseCreator::GetInstance()
{
  return csc;
}

void G4EllipseCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4double semi1,semi2;
  
  // L. Broglia
  // G4Placement* place;
  G4Axis2Placement3D* place=0;

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
  semi1 = *Attr->ptr.r;

  Attr = Ent.NextAttribute();	
  semi2 = *Attr->ptr.r;  
}

void G4EllipseCreator::CreateSTEPGeometry(void* G4obj)
{
}
