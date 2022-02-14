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
// G4LogicalBorderSurface Implementation
//
// A Logical Surface class for surfaces defined by the boundary
// of two physical volumes.
//
// Author: John Apostolakis (John.Apostolakis@cern.ch), 26-06-1997
// --------------------------------------------------------------------

#include "G4LogicalBorderSurface.hh"
#include "G4VPhysicalVolume.hh"

G4LogicalBorderSurfaceTable*
G4LogicalBorderSurface::theBorderSurfaceTable = nullptr;

//
// Constructor & destructor
//

G4LogicalBorderSurface::
G4LogicalBorderSurface(const G4String& name,
                             G4VPhysicalVolume* vol1, 
                             G4VPhysicalVolume* vol2,
                             G4SurfaceProperty* surfaceProperty)
  : G4LogicalSurface(name, surfaceProperty),
    Volume1(vol1), Volume2(vol2),
    Index(theBorderSurfaceTable != nullptr ? theBorderSurfaceTable->size() : 0)
{
  if (theBorderSurfaceTable == nullptr)
  {
    theBorderSurfaceTable = new G4LogicalBorderSurfaceTable;
  }

  // Store in the table of Surfaces
  //
  theBorderSurfaceTable->insert(std::make_pair(std::make_pair(vol1,vol2),this));
}

G4LogicalBorderSurface::~G4LogicalBorderSurface()
{
}

//
// Operators
//

G4bool
G4LogicalBorderSurface::operator==(const G4LogicalBorderSurface &right) const
{
  return (this == (G4LogicalBorderSurface *) &right);
}

G4bool
G4LogicalBorderSurface::operator!=(const G4LogicalBorderSurface &right) const
{
  return (this != (G4LogicalBorderSurface *) &right);
}

//
// Methods
//

const G4LogicalBorderSurfaceTable* G4LogicalBorderSurface::GetSurfaceTable()
{
  if (theBorderSurfaceTable == nullptr)
  {
    theBorderSurfaceTable = new G4LogicalBorderSurfaceTable;
  }
  return theBorderSurfaceTable;
}

std::size_t G4LogicalBorderSurface::GetNumberOfBorderSurfaces()
{
  if (theBorderSurfaceTable != nullptr)
  {
    return theBorderSurfaceTable->size();
  }
  return 0;
}

G4LogicalBorderSurface*
G4LogicalBorderSurface::GetSurface(const G4VPhysicalVolume* vol1,
                                   const G4VPhysicalVolume* vol2)
{
  if (theBorderSurfaceTable != nullptr)
  {
    auto pos = theBorderSurfaceTable->find(std::make_pair(vol1,vol2));
    if(pos != theBorderSurfaceTable->cend()) return pos->second;
  }
  return nullptr;
}

// Dump info for known surfaces
//
void G4LogicalBorderSurface::DumpInfo()
{
  G4cout << "***** Surface Table : Nb of Surfaces = "
         << GetNumberOfBorderSurfaces() << " *****" << G4endl;

  if (theBorderSurfaceTable != nullptr)
  {
    for(auto pos = theBorderSurfaceTable->cbegin();
        pos != theBorderSurfaceTable->cend(); ++pos)
    {
      G4LogicalBorderSurface* pSurf = pos->second;
      G4cout << pSurf->GetName() << " : " << G4endl
             << " Border of volumes "
             << pSurf->GetVolume1()->GetName() << " and " 
             << pSurf->GetVolume2()->GetName() << G4endl;
    }
  }
  G4cout << G4endl;
}

void G4LogicalBorderSurface::CleanSurfaceTable()
{
  if (theBorderSurfaceTable != nullptr)
  {
    for(auto pos = theBorderSurfaceTable->cbegin();
        pos != theBorderSurfaceTable->cend(); ++pos)
    {
      if (pos->second)  { delete pos->second; }
    }
    theBorderSurfaceTable->clear();
  }
  return;
}
