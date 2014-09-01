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
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// 
// G4UTrd
//
// Class description:
//
//   Wrapper class for UTrd to make use of UTrd from USolids module.

// History:
// 13.09.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------
#ifndef G4UTRD_HH
#define G4UTRD_HH

#include "G4USolid.hh"
#include "UTrd.hh"
#include "G4Polyhedron.hh"

class G4UTrd : public G4USolid 
{
  public:  // with description

    G4UTrd(const G4String& pName,
                 G4double pdx1, G4double pdx2,
                 G4double pdy1, G4double pdy2,
                 G4double pdz);
      // Constructs a trapezoid with name, and half lengths

   ~G4UTrd();

    void ComputeDimensions(      G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4VSolid* Clone() const;

    inline UTrd* GetShape() const;

    inline G4double GetXHalfLength1() const;
    inline G4double GetXHalfLength2() const;
    inline G4double GetYHalfLength1() const;
    inline G4double GetYHalfLength2() const;
    inline G4double GetZHalfLength()  const;

    inline void SetXHalfLength1(G4double val);
    inline void SetXHalfLength2(G4double val);
    inline void SetYHalfLength1(G4double val);
    inline void SetYHalfLength2(G4double val);
    inline void SetZHalfLength(G4double val);

    inline void SetAllParameters(G4double pdx1, G4double pdx2,
                                 G4double pdy1, G4double pdy2,
                                 G4double pdz);

  public:  // without description

    G4UTrd(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UTrd(const G4UTrd& rhs);
    G4UTrd& operator=(const G4UTrd& rhs); 
      // Copy constructor and assignment operator.
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline UTrd* G4UTrd::GetShape() const
{
  return (UTrd*) fShape;
}

inline G4double G4UTrd::GetXHalfLength1() const
{
  return GetShape()->GetXHalfLength1();
}
inline G4double G4UTrd::GetXHalfLength2() const
{
  return GetShape()->GetXHalfLength2();
}
inline G4double G4UTrd::GetYHalfLength1() const
{
  return GetShape()->GetYHalfLength1();
}
inline G4double G4UTrd::GetYHalfLength2() const
{
  return GetShape()->GetYHalfLength2();
}
inline G4double G4UTrd::GetZHalfLength()  const
{
  return GetShape()->GetZHalfLength();
}

inline void G4UTrd::SetXHalfLength1(G4double val)
{
  GetShape()->SetXHalfLength1(val);
  fRebuildPolyhedron = true;
}
inline void G4UTrd::SetXHalfLength2(G4double val)
{
  GetShape()->SetXHalfLength2(val);
  fRebuildPolyhedron = true;
}
inline void G4UTrd::SetYHalfLength1(G4double val)
{
  GetShape()->SetYHalfLength1(val);
  fRebuildPolyhedron = true;
}
inline void G4UTrd::SetYHalfLength2(G4double val)
{
  GetShape()->SetYHalfLength2(val);
  fRebuildPolyhedron = true;
}
inline void G4UTrd::SetZHalfLength(G4double val)
{
  GetShape()->SetZHalfLength(val);
  fRebuildPolyhedron = true;
}
inline void G4UTrd::SetAllParameters(G4double pdx1, G4double pdx2,
                                     G4double pdy1, G4double pdy2,
                                     G4double pdz)
{
  GetShape()->SetAllParameters(pdx1, pdx2, pdy1, pdy2, pdz);
  fRebuildPolyhedron = true;
}

#endif
