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
// * This  code  implementation is the  intellectual property  of the *
// * Vanderbilt University Free Electron Laser Center                 *
// * Vanderbilt University, Nashville, TN, USA                        *
// * Development supported by:                                        *
// * United States MFEL program  under grant FA9550-04-1-0045         *
// * and NASA under contract number NNG04CT05P.                       *
// * Written by Marcus H. Mendenhall and Robert A. Weller.            *
// *                                                                  *
// * Contributed to the Geant4 Core, January, 2005.                   *
// *                                                                  *
// ********************************************************************
//
//
// $Id:$
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4UTet
//
// Class description:
//
//   Wrapper class for UTet to make use of UTet from USolids module.

// History:
// 1.11.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------
#ifndef G4UTET_HH
#define G4UTET_HH

#include "G4USolid.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "UTet.hh"

#include "G4Polyhedron.hh"

class G4UTet : public G4USolid
{

  public:  // with description

    G4UTet(const G4String& pName, 
                 G4ThreeVector anchor,
                 G4ThreeVector p2,
                 G4ThreeVector p3,
                 G4ThreeVector p4, 
                 G4bool *degeneracyFlag=0);

   ~G4UTet();

    inline UTet* GetShape() const;

    inline G4GeometryType GetEntityType() const;

  public:   // without description

    G4UTet(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UTet(const G4UTet& rhs);
    G4UTet& operator=(const G4UTet& rhs); 
      // Copy constructor and assignment operator.

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const;

    G4Polyhedron* CreatePolyhedron() const;

    std::vector<G4ThreeVector> GetVertices() const;
      // Return the four vertices of the shape.
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline UTet* G4UTet::GetShape() const
{
  return (UTet*) fShape;
}

inline G4GeometryType G4UTet::GetEntityType() const
{
  return "G4Tet";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
