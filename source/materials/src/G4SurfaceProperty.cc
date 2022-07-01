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

G4SurfaceProperty::G4SurfaceProperty(const G4String& name, G4SurfaceType type)
  : theName(name)
  , theType(type)
{
  theSurfacePropertyTable.push_back(this);
}

G4SurfaceProperty::G4SurfaceProperty()
  : theName("Dielectric")
  , theType(dielectric_metal)
{
  theSurfacePropertyTable.push_back(this);
}

const G4SurfacePropertyTable* G4SurfaceProperty::GetSurfacePropertyTable()
{
  return &theSurfacePropertyTable;
}

size_t G4SurfaceProperty::GetNumberOfSurfaceProperties()
{
  return theSurfacePropertyTable.size();
}

// Dump info for known surface properties
void G4SurfaceProperty::DumpTableInfo()
{
  G4cout << "***** Surface Property Table : Nb of Surface Properties = "
         << GetNumberOfSurfaceProperties() << " *****" << G4endl;

  for(auto pSurfaceProperty : theSurfacePropertyTable)
  {
    G4cout << pSurfaceProperty->GetName() << " : " << G4endl
           << "  Surface Property type   = " << pSurfaceProperty->GetType()
           << G4endl;
  }
  G4cout << G4endl;
}

void G4SurfaceProperty::CleanSurfacePropertyTable()
{
  DumpTableInfo();
  G4SurfacePropertyTable::iterator pos;
  for(pos = theSurfacePropertyTable.begin();
      pos != theSurfacePropertyTable.end(); pos++)
  {
    if(*pos != nullptr)
    {
      delete *pos;
    }
  }
  theSurfacePropertyTable.clear();
  DumpTableInfo();
}
