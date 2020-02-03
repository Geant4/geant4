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
// G4LogicalSkinSurface
//
// Class description:
//
// A Logical Surface class for the surface surrounding a single logical
// volume.

// Author: John Apostolakis (John.Apostolakis@cern.ch), 16-06-1997
//
// ----------------------------------------------------------------------
#ifndef G4LogicalSkinSurface_h
#define G4LogicalSkinSurface_h 1

#include <vector>

#include "G4LogicalSurface.hh"

class G4LogicalVolume;
class G4LogicalSkinSurface;

using G4LogicalSkinSurfaceTable = std::vector<G4LogicalSkinSurface*>;

class G4LogicalSkinSurface : public G4LogicalSurface 
{

  public:  // with description

    G4LogicalSkinSurface( const G4String& name,
                                G4LogicalVolume* vol,
                                G4SurfaceProperty* surfaceProperty );
    ~G4LogicalSkinSurface();
      // Constructor and destructor.

    G4LogicalSkinSurface(const G4LogicalSkinSurface&) = delete;
    G4LogicalSkinSurface& operator=(const G4LogicalSkinSurface&) = delete;
      // Assignment and copying not allowed.

    G4bool operator==(const G4LogicalSkinSurface &right) const;
    G4bool operator!=(const G4LogicalSkinSurface &right) const;
      // Operators.

    static G4LogicalSkinSurface* GetSurface(const G4LogicalVolume* vol);
    inline const G4LogicalVolume* GetLogicalVolume() const;
    inline void  SetLogicalVolume(G4LogicalVolume* vol);
      // Accessors.

    static void CleanSurfaceTable();
    static const G4LogicalSkinSurfaceTable* GetSurfaceTable();
    static size_t GetNumberOfSkinSurfaces();
    static void DumpInfo(); // const 
      // To handle with the table of surfaces.

  private:

    G4LogicalVolume* LogVolume;
      // Logical Volume pointer on side 1.

    static G4LogicalSkinSurfaceTable *theSkinSurfaceTable;
      // The static Table of SkinSurfaces.

};

// ********************************************************************
// Inline methods
// ********************************************************************

#include "G4LogicalSkinSurface.icc"

#endif /* G4LogicalSkinSurface_h */

