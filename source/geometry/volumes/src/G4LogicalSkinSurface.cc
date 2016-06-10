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
// $Id: G4LogicalSkinSurface.cc 87615 2014-12-12 15:23:11Z gcosmo $
//
// --------------------------------------------------------------------
// G4LogicalSkinSurface Implementation
// --------------------------------------------------------------------
//
// A Logical Surface class for the surface surrounding a single
// logical volume.
//
// Created:     1997-06-26
// Author:      John Apostolakis (John.Apostolakis@cern.ch)
//
// ----------------------------------------------------------------------

#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"

G4LogicalSkinSurfaceTable *G4LogicalSkinSurface::theSkinSurfaceTable = 0;

//
// Constructors & destructor
//

G4LogicalSkinSurface::G4LogicalSkinSurface(const G4String&   name,
                                           G4LogicalVolume*  logicalVolume,
                                           G4SurfaceProperty* surfaceProperty)
  : G4LogicalSurface(name, surfaceProperty),
    LogVolume(logicalVolume)
{
  if (!theSkinSurfaceTable)
  {
    theSkinSurfaceTable = new G4LogicalSkinSurfaceTable;
  }
  // Store in the table of Surfaces
  //
  theSkinSurfaceTable->push_back(this);
}

G4LogicalSkinSurface::G4LogicalSkinSurface(const G4LogicalSkinSurface& right)
  : G4LogicalSurface(right.GetName(), right.GetSurfaceProperty())
{
  if (!theSkinSurfaceTable)
  {
    theSkinSurfaceTable = new G4LogicalSkinSurfaceTable;
  }
  SetTransitionRadiationSurface(right.GetTransitionRadiationSurface());
  LogVolume = right.LogVolume;
  (*theSkinSurfaceTable) = (*right.theSkinSurfaceTable);
}

G4LogicalSkinSurface::~G4LogicalSkinSurface()
{
//  delete theSkinSurfaceTable; theSkinSurfaceTable = 0;
}

//
// Operators
//

G4LogicalSkinSurface&
G4LogicalSkinSurface::operator=(const G4LogicalSkinSurface& right)
{
  if (&right == this) return *this;

  SetSurfaceProperty(right.GetSurfaceProperty());
  SetName(right.GetName());
  SetTransitionRadiationSurface(right.GetTransitionRadiationSurface());
  LogVolume = right.LogVolume;
  (*theSkinSurfaceTable) = (*right.theSkinSurfaceTable);

  return *this;
}

G4int
G4LogicalSkinSurface::operator==(const G4LogicalSkinSurface& right) const
{
  return (this == (G4LogicalSkinSurface *) &right);
}

G4int
G4LogicalSkinSurface::operator!=(const G4LogicalSkinSurface& right) const
{
  return (this != (G4LogicalSkinSurface *) &right);
}

//
// Methods
//

const G4LogicalSkinSurfaceTable* G4LogicalSkinSurface::GetSurfaceTable()
{
  if (!theSkinSurfaceTable)
  {
    theSkinSurfaceTable = new G4LogicalSkinSurfaceTable;
  }
  return theSkinSurfaceTable;
}

size_t G4LogicalSkinSurface::GetNumberOfSkinSurfaces()
{
  if (theSkinSurfaceTable)
  {
    return theSkinSurfaceTable->size();
  }
  return 0;
}

G4LogicalSkinSurface*
G4LogicalSkinSurface::GetSurface(const G4LogicalVolume* vol)
{
  if (theSkinSurfaceTable)
  {
    for (size_t i=0; i<theSkinSurfaceTable->size(); i++)
    {
      if((*theSkinSurfaceTable)[i]->GetLogicalVolume() == vol)
        return (*theSkinSurfaceTable)[i];
    }
  }
  return 0;
}

// Dump info for known surfaces
//
void G4LogicalSkinSurface::DumpInfo() 
{
  G4cout << "***** Skin Surface Table : Nb of Surfaces = "
         << GetNumberOfSkinSurfaces() << " *****" << G4endl;

  if (theSkinSurfaceTable)
  {
    for (size_t i=0; i<theSkinSurfaceTable->size(); i++)
    {
      G4LogicalSkinSurface* pSkinSurface = (*theSkinSurfaceTable)[i];
      G4cout << pSkinSurface->GetName() << " : " << G4endl
             << " Skin of logical volume "
             << pSkinSurface->GetLogicalVolume()->GetName()
             << G4endl;
    }
  }
  G4cout << G4endl;
}

void G4LogicalSkinSurface::CleanSurfaceTable()
{
  if (theSkinSurfaceTable)
  {
    G4LogicalSkinSurfaceTable::iterator pos;
    for(pos=theSkinSurfaceTable->begin();
        pos!=theSkinSurfaceTable->end(); pos++)
    {
      if (*pos) { delete *pos; }
    }
    theSkinSurfaceTable->clear();
  }
  return;
}
