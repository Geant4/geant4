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
// $Id: G4GeometricRepresentationContextCreator.cc,v 1.5 2001-07-11 10:00:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4GeometricRepresentationContextCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4GeometricRepresentationContextCreator.hh"
#include "G4GeometryTable.hh"

G4GeometricRepresentationContextCreator
  G4GeometricRepresentationContextCreator::csc;

G4GeometricRepresentationContextCreator::
  G4GeometricRepresentationContextCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4GeometricRepresentationContextCreator::
  ~G4GeometricRepresentationContextCreator() {}

void G4GeometricRepresentationContextCreator::CreateG4Geometry(STEPentity& Ent)
{
  // Made by L. Broglia
  STEPattribute *Attr;

  G4String attrName("coordinate_space_dimension");
  Attr = GetNamedAttribute(attrName, Ent);

  // Get coordinate space dimension (ptr.r --> REAL_TYPE)
  // SdaiReal Tmpdim = *Attr->ptr.r;

  // void *tmp = G4GeometryTable::CreateObject(*TmpEnt);
  // place = (G4Axis2Placement3D*)tmp...;

}

void G4GeometricRepresentationContextCreator::CreateSTEPGeometry(void* G4obj)
{
}
