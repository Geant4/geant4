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
// $Id: G4LogicalBorderSurface.cc 87615 2014-12-12 15:23:11Z gcosmo $
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
// ----------------------------------------------------------------------

#include "G4LogicalBorderSurface.hh"
#include "G4VPhysicalVolume.hh"

G4LogicalBorderSurfaceTable *G4LogicalBorderSurface::theBorderSurfaceTable = 0;

//
// Constructor & destructor
//

G4LogicalBorderSurface::
G4LogicalBorderSurface(const G4String& name,
                             G4VPhysicalVolume* vol1, 
                             G4VPhysicalVolume* vol2,
                             G4SurfaceProperty* surfaceProperty)
  : G4LogicalSurface(name, surfaceProperty),
    Volume1(vol1), Volume2(vol2)
{
  if (!theBorderSurfaceTable)
  {
    theBorderSurfaceTable = new G4LogicalBorderSurfaceTable;
  }

  // Store in the table of Surfaces
  //
  theBorderSurfaceTable->push_back(this);
}

G4LogicalBorderSurface::
G4LogicalBorderSurface(const G4LogicalBorderSurface& right)
  : G4LogicalSurface(right.GetName(), right.GetSurfaceProperty())
{
  if (!theBorderSurfaceTable)
  {
    theBorderSurfaceTable = new G4LogicalBorderSurfaceTable;
  }
  SetTransitionRadiationSurface(right.GetTransitionRadiationSurface());
  Volume1 = right.Volume1;
  Volume2 = right.Volume2;
  (*theBorderSurfaceTable) = (*right.theBorderSurfaceTable);
}

G4LogicalBorderSurface::~G4LogicalBorderSurface()
{
//  delete theBorderSurfaceTable; theBorderSurfaceTable=0;
}

//
// Operators
//

G4LogicalBorderSurface&
G4LogicalBorderSurface::operator=(const G4LogicalBorderSurface &right)
{
  if (&right == this) return *this;

  SetSurfaceProperty(right.GetSurfaceProperty());
  SetName(right.GetName());
  SetTransitionRadiationSurface(right.GetTransitionRadiationSurface());
  Volume1 = right.Volume1;
  Volume2 = right.Volume2;
  (*theBorderSurfaceTable) = (*right.theBorderSurfaceTable);

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

const G4LogicalBorderSurfaceTable* G4LogicalBorderSurface::GetSurfaceTable()
{
  if (!theBorderSurfaceTable)
  {
    theBorderSurfaceTable = new G4LogicalBorderSurfaceTable;
  }
  return theBorderSurfaceTable;
}

size_t G4LogicalBorderSurface::GetNumberOfBorderSurfaces()
{
  if (theBorderSurfaceTable)
  {
    return theBorderSurfaceTable->size();
  }
  return 0;
}

G4LogicalBorderSurface*
G4LogicalBorderSurface::GetSurface(const G4VPhysicalVolume* vol1,
                                   const G4VPhysicalVolume* vol2)
{
  if (theBorderSurfaceTable)
  {
    for (size_t i=0; i<theBorderSurfaceTable->size(); i++)
    {
      if( ((*theBorderSurfaceTable)[i]->GetVolume1() == vol1) &&
          ((*theBorderSurfaceTable)[i]->GetVolume2() == vol2) )
        return (*theBorderSurfaceTable)[i];
    }
  }
  return 0;
}

// Dump info for known surfaces
//
void G4LogicalBorderSurface::DumpInfo()
{
  G4cout << "***** Surface Table : Nb of Surfaces = "
         << GetNumberOfBorderSurfaces() << " *****" << G4endl;

  if (theBorderSurfaceTable)
  {
    for (size_t i=0; i<theBorderSurfaceTable->size(); i++)
    {
      G4LogicalBorderSurface* pBorderSurface = (*theBorderSurfaceTable)[i];
      G4cout << pBorderSurface->GetName() << " : " << G4endl
             << " Border of volumes "
             << pBorderSurface->GetVolume1()->GetName() << " and " 
             << pBorderSurface->GetVolume2()->GetName()
             << G4endl;
    }
  }
  G4cout << G4endl;
}

void G4LogicalBorderSurface::CleanSurfaceTable()
{
  if (theBorderSurfaceTable)
  {
    G4LogicalBorderSurfaceTable::iterator pos;
    for(pos=theBorderSurfaceTable->begin();
        pos!=theBorderSurfaceTable->end(); pos++)
    {
      if (*pos)  { delete *pos; }
    }
    theBorderSurfaceTable->clear();
  }
  return;
}
