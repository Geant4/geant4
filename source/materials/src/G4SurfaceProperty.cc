//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4SurfaceProperty.cc,v 1.1 2007/04/25 16:19:26 gum Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// 
////////////////////////////////////////////////////////////////////////
// G4SurfaceProperty Implementation
////////////////////////////////////////////////////////////////////////
//
// Class Description:
//
// A base class describing a surface property.
// Derived classes are G4Opticalsurface, G4Firovsurface, etc.      
//
// File:        G4SurfaceProperty.cc
// Version:     1.0
// Created:     16-11-2006
// Author:      Peter Gumplinger
//
////////////////////////////////////////////////////////////////////////

#include "globals.hh"

#include "G4SurfaceProperty.hh"

G4SurfacePropertyTable G4SurfaceProperty::theSurfacePropertyTable;

//
// Methods
//

const G4SurfacePropertyTable* G4SurfaceProperty::GetSurfacePropertyTable()
{
  return &theSurfacePropertyTable;
}

size_t G4SurfaceProperty::GetNumberOfSurfaceProperties()
{
  return theSurfacePropertyTable.size();
}

// Dump info for known surface properties
//
void G4SurfaceProperty::DumpInfo()
{
  G4cout << "***** Surface Property Table : Nb of Surface Properties = "
         << GetNumberOfSurfaceProperties() << " *****" << G4endl;

  for (size_t i=0; i<theSurfacePropertyTable.size(); i++)
  {
    G4SurfaceProperty* pSurfaceProperty = theSurfacePropertyTable[i];
    G4cout << pSurfaceProperty->GetName() << " : " << G4endl
           << "  Surface Property type   = " 
           << pSurfaceProperty->GetType()
           << G4endl;
  }
  G4cout << G4endl;
}

void G4SurfaceProperty::CleanSurfacePropertyTable()
{
  DumpInfo();
  G4SurfacePropertyTable::iterator pos;
  for(pos=theSurfacePropertyTable.begin();
      pos!=theSurfacePropertyTable.end(); pos++)
  {
    if (*pos) delete *pos;
  }
  theSurfacePropertyTable.clear();
  DumpInfo();
}
