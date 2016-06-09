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
// $Id: G4LogicalSkinSurface.cc,v 1.12 2004/05/19 08:14:42 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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
#include "G4ios.hh"
#include "G4LogicalVolume.hh"

G4LogicalSkinSurfaceTable G4LogicalSkinSurface::theSurfaceTable;

//
// Constructors & destructor
//

G4LogicalSkinSurface::G4LogicalSkinSurface(const G4String&   name,
                                           G4LogicalVolume*  logicalVolume,
                                           G4SurfaceProperty* surfaceProperty)
  : G4LogicalSurface(name, surfaceProperty),
    LogVolume(logicalVolume)
{
  // Store in the table of Surfaces
  //
  theSurfaceTable.push_back(this);
}

G4LogicalSkinSurface::G4LogicalSkinSurface(const G4LogicalSkinSurface& right)
  : G4LogicalSurface(right.GetName(), right.GetSurfaceProperty())
{
  SetTransitionRadiationSurface(right.GetTransitionRadiationSurface());
  LogVolume = right.LogVolume;
  theSurfaceTable = right.theSurfaceTable;
}

G4LogicalSkinSurface::~G4LogicalSkinSurface()
{
}

//
// Operators
//

const G4LogicalSkinSurface&
G4LogicalSkinSurface::operator=(const G4LogicalSkinSurface& right)
{
  if (&right == this) return *this;
  if (&right)
  {
    SetSurfaceProperty(right.GetSurfaceProperty());
    SetName(right.GetName());
    SetTransitionRadiationSurface(right.GetTransitionRadiationSurface());
    LogVolume = right.LogVolume;
    theSurfaceTable = right.theSurfaceTable;
  }
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

size_t G4LogicalSkinSurface::GetNumberOfSkinSurfaces()
{
  return theSurfaceTable.size();
}

G4LogicalSkinSurface*
G4LogicalSkinSurface::GetSurface(const G4LogicalVolume* vol)
{
  for (size_t i=0; i<theSurfaceTable.size(); i++)
  {
    if(theSurfaceTable[i]->GetLogicalVolume() == vol)
      return theSurfaceTable[i];
  }
  return NULL;
}

// Dump info for known surfaces
//
void G4LogicalSkinSurface::DumpInfo() 
{
    G4cout << "***** Surface Table : Nb of Surfaces = "
           << GetNumberOfSkinSurfaces() << " *****" << G4endl;

    for (size_t i=0; i<theSurfaceTable.size(); i++)
    {
      G4LogicalSkinSurface* pSkinSurface = theSurfaceTable[i];
      G4cout << theSurfaceTable[i]->GetName() << " : " << G4endl
             << " Skin of logical volume "
             << pSkinSurface->GetLogicalVolume()->GetName()
             << G4endl;
    }
    G4cout << G4endl;
}

void G4LogicalSkinSurface::CleanSurfaceTable()
{  
  std::vector<G4LogicalSkinSurface*>::iterator pos;
  for(pos=theSurfaceTable.begin(); pos!=theSurfaceTable.end(); pos++)
  {
    if (*pos) delete *pos;
  }
  theSurfaceTable.clear();
}
