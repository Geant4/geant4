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
// $Id: G4LogicalBorderSurface.cc,v 1.8 2002-08-06 08:23:38 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
// G4LogicalBorderSurface Implementation
// --------------------------------------------------------------------
//
// A Logical Surface class for surfaces defined by the boundary
// of two physical volumes.
//
// Created:     1997-06-26
// Author:      John Apostolakis (John.Apostolakis@cern.ch)
//
// ********************************************************************

#include "G4LogicalBorderSurface.hh"

G4LogicalBorderSurfaceTable G4LogicalBorderSurface::theBorderSurfaceTable;

//
// Constructor & destructor
//

G4LogicalBorderSurface::G4LogicalBorderSurface(const G4String& name,
                                               G4VPhysicalVolume* vol1, 
                                               G4VPhysicalVolume* vol2,
                                               G4OpticalSurface* opticsSurface)
  : G4LogicalSurface(name, opticsSurface),
    Volume1(vol1), Volume2(vol2)
{
  // Store in the table of Surfaces
  //
  theBorderSurfaceTable.push_back(this);
}

G4LogicalBorderSurface::
G4LogicalBorderSurface(const G4LogicalBorderSurface& right)
  : G4LogicalSurface(right.GetName(), right.GetOpticalSurface())
{
  SetTransitionRadiationSurface(right.GetTransitionRadiationSurface());
  Volume1 = right.Volume1;
  Volume2 = right.Volume2;
  theBorderSurfaceTable = right.theBorderSurfaceTable;
}

G4LogicalBorderSurface::~G4LogicalBorderSurface()
{
}

//
// Operators
//

const G4LogicalBorderSurface&
G4LogicalBorderSurface::operator=(const G4LogicalBorderSurface &right)
{
  if (&right == this) return *this;
  if (&right)
  {
    SetOpticalSurface(right.GetOpticalSurface());
    SetName(right.GetName());
    SetTransitionRadiationSurface(right.GetTransitionRadiationSurface());
    Volume1 = right.Volume1;
    Volume2 = right.Volume2;
    theBorderSurfaceTable = right.theBorderSurfaceTable;
  }
  return *this;
}

G4int
G4LogicalBorderSurface::operator==(const G4LogicalBorderSurface &right) const
{
  return (this == (G4LogicalBorderSurface *) &right);
}

G4int
G4LogicalBorderSurface::operator!=(const G4LogicalBorderSurface &right) const
{
  return (this != (G4LogicalBorderSurface *) &right);
}

//
// Methods
//

const G4std::vector<G4LogicalBorderSurface*> *
G4LogicalBorderSurface::GetSurfaceTable()
{
  return &theBorderSurfaceTable;
}

size_t G4LogicalBorderSurface::GetNumberOfBorderSurfaces()
{
  return theBorderSurfaceTable.size();
}

G4LogicalBorderSurface*
G4LogicalBorderSurface::GetSurface(const G4VPhysicalVolume* vol1,
                                   const G4VPhysicalVolume* vol2)
{
  for (size_t i=0; i<theBorderSurfaceTable.size(); i++)
  {
    if( (theBorderSurfaceTable[i]->GetVolume1() == vol1) &&
        (theBorderSurfaceTable[i]->GetVolume2() == vol2) )
      return theBorderSurfaceTable[i];
  }
  return NULL;
}

// Dump info for known surfaces
//
void G4LogicalBorderSurface::DumpInfo()
{
  G4cout << "***** Surface Table : Nb of Surfaces = "
         << GetNumberOfBorderSurfaces() << " *****" << G4endl;

  for (size_t i=0; i<theBorderSurfaceTable.size(); i++)
  {
    G4cout << theBorderSurfaceTable[i]->GetName() << " : " << G4endl
           << "  Surface type   = " << theBorderSurfaceTable[i]->GetName()
           << G4endl;
  }
  G4cout << G4endl;
}
