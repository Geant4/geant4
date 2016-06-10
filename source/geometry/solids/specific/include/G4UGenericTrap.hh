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
// G4UGenericTrap
//
// Class description:
//
//   Wrapper class for G4UGenericTrap to make use of it from USolids module.

// History:
// 30.10.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------
#ifndef G4UGENERICTRAP_hh
#define G4UGENERICTRAP_hh

#include "G4USolid.hh"
#include "UGenericTrap.hh"
#include "G4TwoVector.hh"

class G4Polyhedron;

class G4UGenericTrap : public G4USolid 
{
  public:  // with description

    G4UGenericTrap(const G4String& name, G4double halfZ,
                   const std::vector<G4TwoVector>& vertices);

   ~G4UGenericTrap();

    inline UGenericTrap* GetShape() const;

    inline G4double    GetZHalfLength() const;
    inline G4int       GetNofVertices() const;
    inline G4TwoVector GetVertex(G4int index) const;
    inline const std::vector<G4TwoVector>& GetVertices() const;
    inline G4double    GetTwistAngle(G4int index) const;
    inline G4bool      IsTwisted() const;
    inline G4int       GetVisSubdivisions() const;
    inline void        SetVisSubdivisions(G4int subdiv);
    inline void        SetZHalfLength(G4double);

  public:  // without description

    G4UGenericTrap(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UGenericTrap( const G4UGenericTrap& source );
    G4UGenericTrap &operator=(const G4UGenericTrap& source);
      // Copy constructor and assignment operator.

    G4Polyhedron* CreatePolyhedron() const;
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline UGenericTrap* G4UGenericTrap::GetShape() const
{
  return (UGenericTrap*) fShape;
}

inline G4double G4UGenericTrap::GetZHalfLength() const
{
  return GetShape()->GetZHalfLength();
}
inline G4int G4UGenericTrap::GetNofVertices() const
{
  return GetShape()->GetNofVertices();
}
inline G4TwoVector G4UGenericTrap::GetVertex(G4int index) const
{
  UVector2 v = GetShape()->GetVertex(index);
  return G4TwoVector(v.x, v.y);
}
inline const std::vector<G4TwoVector>& G4UGenericTrap::GetVertices() const
{
  std::vector<UVector2> v = GetShape()->GetVertices();
  static std::vector<G4TwoVector> vertices; vertices.clear();
  for (size_t n=0; n<v.size(); ++n)
  {
    vertices.push_back(G4TwoVector(v[n].x,v[n].y));
  }
  return vertices;
}
inline G4double G4UGenericTrap::GetTwistAngle(G4int index) const
{
  return GetShape()->GetTwistAngle(index);
}
inline G4bool G4UGenericTrap::IsTwisted() const
{
  return GetShape()->IsTwisted();
}
inline G4int G4UGenericTrap::GetVisSubdivisions() const
{
  return GetShape()->GetVisSubdivisions();
}

inline void G4UGenericTrap::SetVisSubdivisions(G4int subdiv)
{
  GetShape()->SetVisSubdivisions(subdiv);
}

inline void G4UGenericTrap::SetZHalfLength(G4double halfZ)
{
  GetShape()->SetZHalfLength(halfZ);
}

#endif
