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
// $Id: G4FaceOuterBoundCreator.cc,v 1.5 2003/03/28 15:47:42 gcosmo Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// ----------------------------------------------------------------------
// Class G4FaceOuterBoundCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4FaceOuterBoundCreator.hh"
#include "G4GeometryTable.hh"

G4FaceOuterBoundCreator G4FaceOuterBoundCreator::csc;

G4FaceOuterBoundCreator G4FaceOuterBoundCreator::GetInstance()
{
  return csc;
}

G4FaceOuterBoundCreator::G4FaceOuterBoundCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4FaceOuterBoundCreator::~G4FaceOuterBoundCreator() {}

void G4FaceOuterBoundCreator::CreateSTEPGeometry(void*)
{
}
