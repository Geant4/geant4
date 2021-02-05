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
// G4LogicalBorderSurface
//
// Class description:
//
// A Logical Surface class for surfaces defined by the boundary
// of two physical volumes.

// Author: John Apostolakis (John.Apostolakis@cern.ch), 17-06-1997
// --------------------------------------------------------------------
#ifndef G4LogicalBorderSurface_hh
#define G4LogicalBorderSurface_hh 1

#include <map>

#include "G4LogicalSurface.hh"
#include "G4VPhysicalVolume.hh"

class G4VPhysicalVolume;
class G4LogicalBorderSurface;

using G4LogicalBorderSurfaceTable
      = std::map<std::pair<const G4VPhysicalVolume*,
                           const G4VPhysicalVolume*>, G4LogicalBorderSurface*>;

class G4LogicalBorderSurface : public G4LogicalSurface
{

  public:

    G4LogicalBorderSurface( const G4String& name,
                                  G4VPhysicalVolume* vol1, 
                                  G4VPhysicalVolume* vol2,
                                  G4SurfaceProperty* surfaceProperty );
    ~G4LogicalBorderSurface();
      // Constructor and destructor

    G4LogicalBorderSurface(const G4LogicalBorderSurface&) = delete;
    G4LogicalBorderSurface& operator=(const G4LogicalBorderSurface&) = delete;
      // Copy constructor and assignment operator not allowed.

    G4bool operator==( const G4LogicalBorderSurface& right ) const;
    G4bool operator!=( const G4LogicalBorderSurface& right ) const;
      // Operators

    static G4LogicalBorderSurface* GetSurface( const G4VPhysicalVolume* vol1,
                                               const G4VPhysicalVolume* vol2 );
    inline void SetPhysicalVolumes( G4VPhysicalVolume* vol1,
                                    G4VPhysicalVolume* vol2 );
    inline const G4VPhysicalVolume* GetVolume1() const;
    inline const G4VPhysicalVolume* GetVolume2() const;
    inline std::size_t GetIndex() const;
      // Generic accessors

    inline void SetVolume1( G4VPhysicalVolume* vol1 );
    inline void SetVolume2( G4VPhysicalVolume* vol2 );
      // To use with care!

    static void CleanSurfaceTable();
    static const G4LogicalBorderSurfaceTable* GetSurfaceTable();
    static std::size_t GetNumberOfBorderSurfaces();
    static void DumpInfo(); 
      // To handle the table of surfaces

  private:

    G4VPhysicalVolume* Volume1;  // Physical Volume pointer on side 1
    G4VPhysicalVolume* Volume2;  // Physical Volume pointer on side 2

    std::size_t Index;           // Creation order index

    static G4LogicalBorderSurfaceTable *theBorderSurfaceTable;
      // The static Table of BorderSurfaces
};

// ********************************************************************
// Inline methods
// ********************************************************************

#include "G4LogicalBorderSurface.icc"

#endif
