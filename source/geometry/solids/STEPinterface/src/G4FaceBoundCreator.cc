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
// $Id: G4FaceBoundCreator.cc,v 1.6 2002-11-21 16:49:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4FaceBoundCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4FaceBoundCreator.hh"
#include "G4GeometryTable.hh"

G4FaceBoundCreator G4FaceBoundCreator::csc;

G4FaceBoundCreator::G4FaceBoundCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4FaceBoundCreator::~G4FaceBoundCreator() {}

G4FaceBoundCreator G4FaceBoundCreator::GetInstance()
{
  return csc;
}

void G4FaceBoundCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4String attrName("bound");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);

  // Get curve
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  if (!tmp)
    G4cerr << "WARNING - G4FaceBoundCreator::CreateG4Geometry" << G4endl
           << "\tUnexpected NULL curve vector (G4CurveVector) !" << G4endl
	   << "\tEntity NOT created." << G4endl;

  // L. Broglia
  // Mistake : the created object tmp is a G4CurveVector
  // G4Curve* crv = (G4Curve*)tmp; 
  G4CurveVector* crv = (G4CurveVector*)tmp; 

  Attr = Ent.NextAttribute();
  //orientation = *Attr->ptr.i;
  //crv->SetSameSense(orientation);
  createdObject = crv;
}

void G4FaceBoundCreator::CreateSTEPGeometry(void* G4obj)
{
}
