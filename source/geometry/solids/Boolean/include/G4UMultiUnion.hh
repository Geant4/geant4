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
// $Id:$
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4UMultiUnion
//
// Class description:
//
//   Wrapper class for G4UMultiUnion to make use of it from USolids module.

// History:
// 30.10.13 G.Cosmo, T.Nikitina, CERN/PH
// --------------------------------------------------------------------
#ifndef G4UMULTIUNION_HH
#define G4UMULTIUNION_HH

#include <CLHEP/Vector/Rotation.h>

#if defined(G4GEOM_USE_USOLIDS)

#include "G4USolid.hh"
#include "UMultiUnion.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "HepPolyhedronProcessor.h"

class G4UMultiUnion : public G4USolid 
{
  public:  // with description

    G4UMultiUnion(const G4String& name);
   ~G4UMultiUnion();

    inline UMultiUnion* GetShape() const;

    void AddNode(G4VSolid& solid, G4Transform3D& trans);
      // Build the multiple union by adding nodes
    G4Transform3D* GetTransformation(G4int index) const;
    G4VSolid* GetSolid(G4int index) const;
    G4int GetNumberOfSolids()const;
    void Voxelize();

  public:  // without description

    G4UMultiUnion(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UMultiUnion( const G4UMultiUnion& source );
    G4UMultiUnion& operator=(const G4UMultiUnion& source);
      // Copy constructor and assignment operator.

    void Extent(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                           G4double& pMin, G4double& pMax) const;

    G4Polyhedron* CreatePolyhedron() const;
      // Called by visualization engine
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline UMultiUnion* G4UMultiUnion::GetShape() const
{
  return (UMultiUnion*) fShape;
}

#endif  // G4GEOM_USE_USOLIDS

#endif
