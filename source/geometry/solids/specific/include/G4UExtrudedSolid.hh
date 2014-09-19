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
// G4UExtrudedSolid
//
// Class description:
//
//   Wrapper class for G4UExtrudedSolid to make use of it from USolids module.

// History:
// 30.10.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------
#ifndef G4UEXTRUDEDSOLID_hh
#define G4UEXTRUDEDSOLID_hh

#include "G4USolid.hh"
#include "UExtrudedSolid.hh"
#include "G4TwoVector.hh"

class G4UExtrudedSolid : public G4USolid 
{
  public:  // without description

    struct ZSection
    {
      ZSection(G4double z, G4TwoVector offset, G4double scale)
        : fZ(z), fOffset(offset), fScale(scale) {}
      ZSection(const UExtrudedSolid::ZSection& zs)
        : fZ(zs.fZ), fOffset(G4TwoVector(zs.fOffset.x,zs.fOffset.y)),
          fScale(zs.fScale) {}

      G4double    fZ;
      G4TwoVector fOffset;
      G4double    fScale;
    };

  public:  // with description

    G4UExtrudedSolid(const G4String&        pName,
                   std::vector<G4TwoVector> polygon,
                   std::vector<ZSection>    zsections);
    // General constructor

    G4UExtrudedSolid(const G4String&        pName,
                   std::vector<G4TwoVector> polygon,
                   G4double                 halfZ,
                   G4TwoVector off1, G4double scale1,
                   G4TwoVector off2, G4double scale2);
    // Special constructor for solid with 2 z-sections

   ~G4UExtrudedSolid();

    inline UExtrudedSolid* GetShape() const;

    inline G4int GetNofVertices() const;
    inline G4TwoVector GetVertex(G4int index) const;
    inline std::vector<G4TwoVector> GetPolygon() const;
    inline G4int GetNofZSections() const;
    inline ZSection GetZSection(G4int index) const;
    inline std::vector<ZSection> GetZSections() const;
  
  public:  // without description

    G4UExtrudedSolid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UExtrudedSolid( const G4UExtrudedSolid& source );
    G4UExtrudedSolid &operator=(const G4UExtrudedSolid& source);
      // Copy constructor and assignment operator.
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline UExtrudedSolid* G4UExtrudedSolid::GetShape() const
{
  return (UExtrudedSolid*) fShape;
}

inline G4int G4UExtrudedSolid::GetNofVertices() const
{
  return GetShape()->GetNofVertices();
}
inline G4TwoVector G4UExtrudedSolid::GetVertex(G4int i) const
{
  UVector2 v = GetShape()->GetVertex(i);
  return G4TwoVector(v.x, v.y);
}
inline std::vector<G4TwoVector> G4UExtrudedSolid::GetPolygon() const
{
  std::vector<UVector2> pol = GetShape()->GetPolygon();
  std::vector<G4TwoVector> v;
  for (unsigned int i=0; i<pol.size(); ++i)
  {
    v.push_back(G4TwoVector(pol[i].x, pol[i].y));
  }
  return v;
}
inline G4int G4UExtrudedSolid::GetNofZSections() const
{
  return GetShape()->GetNofZSections();
}
inline G4UExtrudedSolid::ZSection G4UExtrudedSolid::GetZSection(G4int i) const
{
  return ZSection(GetShape()->GetZSection(i));
}
inline std::vector<G4UExtrudedSolid::ZSection> G4UExtrudedSolid::GetZSections() const
{
  std::vector<UExtrudedSolid::ZSection> sv = GetShape()->GetZSections();
  std::vector<G4UExtrudedSolid::ZSection> vec;
  for (unsigned int i=0; i<sv.size(); ++i)
  {
    vec.push_back(ZSection(sv[i]));
  }
  return vec;
}
  
#endif
