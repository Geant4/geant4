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
// Class G4LogicalSkinSurface Implementation
//
// Author: John Apostolakis (CERN), 16 June 1997
// --------------------------------------------------------------------

#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4GeometryManager.hh"

G4LogicalSkinSurfaceTable *G4LogicalSkinSurface::theSkinSurfaceTable = nullptr;

// --------------------------------------------------------------------
// Constructor
//
G4LogicalSkinSurface::G4LogicalSkinSurface(const G4String&  name,
                                           G4LogicalVolume* logicalVolume,
                                           G4SurfaceProperty* surfaceProperty)
  : G4LogicalSurface(name, surfaceProperty),
    LogVolume(logicalVolume)
{
  if (theSkinSurfaceTable == nullptr)
  {
    theSkinSurfaceTable = new G4LogicalSkinSurfaceTable;
  }
  // Store in the table of Surfaces
  //
  theSkinSurfaceTable->insert(std::make_pair(logicalVolume, this));
}

// --------------------------------------------------------------------
// Default destructor
//
G4LogicalSkinSurface::~G4LogicalSkinSurface() = default;

// --------------------------------------------------------------------
G4bool
G4LogicalSkinSurface::operator==(const G4LogicalSkinSurface& right) const
{
  return (this == (G4LogicalSkinSurface *) &right);
}

// --------------------------------------------------------------------
G4bool
G4LogicalSkinSurface::operator!=(const G4LogicalSkinSurface& right) const
{
  return (this != (G4LogicalSkinSurface *) &right);
}

// --------------------------------------------------------------------
const G4LogicalSkinSurfaceTable* G4LogicalSkinSurface::GetSurfaceTable()
{
  if (theSkinSurfaceTable == nullptr)
  {
    theSkinSurfaceTable = new G4LogicalSkinSurfaceTable;
  }
  return theSkinSurfaceTable;
}

// --------------------------------------------------------------------
size_t G4LogicalSkinSurface::GetNumberOfSkinSurfaces()
{
  if (theSkinSurfaceTable != nullptr)
  {
    return theSkinSurfaceTable->size();
  }
  return 0;
}

// --------------------------------------------------------------------
G4LogicalSkinSurface*
G4LogicalSkinSurface::GetSurface(const G4LogicalVolume* vol)
{
  if (theSkinSurfaceTable != nullptr)
  {
    auto pos = theSkinSurfaceTable->find(vol);
    if(pos != theSkinSurfaceTable->cend()) { return pos->second; }
  }
  return nullptr;
}

// --------------------------------------------------------------------
// Dump info for known surfaces
//
void G4LogicalSkinSurface::DumpInfo() 
{
  G4cout << "***** Skin Surface Table : Nb of Surfaces = "
         << GetNumberOfSkinSurfaces() << " *****" << G4endl;

  if (theSkinSurfaceTable != nullptr)
  {
    for(const auto & pos : *theSkinSurfaceTable)
    {
      G4LogicalSkinSurface* pSurf = pos.second;
      G4cout << pSurf->GetName() << " : " << G4endl
             << " Skin of logical volume "
             << pSurf->GetLogicalVolume()->GetName()
             << G4endl;
    }
  }
  G4cout << G4endl;
}

// --------------------------------------------------------------------
void G4LogicalSkinSurface::CleanSurfaceTable()
{
  if (theSkinSurfaceTable == nullptr) { return; }

  // Do nothing if geometry is closed
  if (G4GeometryManager::GetInstance()->IsGeometryClosed())
  {
    G4cout << "WARNING - Attempt to clear the skin surface store"
           << " while geometry closed !" << G4endl;
    return;
  }

  for (const auto& pos : *theSkinSurfaceTable)
  {
    delete pos.second;
  }
  theSkinSurfaceTable->clear();
}
