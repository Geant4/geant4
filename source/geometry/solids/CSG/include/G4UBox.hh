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
// G4UBox
//
// Class description:
//
//   Wrapper class for UBox to make use of UBox from USolids module.

// History:
// 13.09.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------
#ifndef G4UBOX_HH
#define G4UBOX_HH

#include "G4USolid.hh"
#include "UBox.hh"
#include "G4Polyhedron.hh"

class G4UBox : public G4USolid 
{
  public:  // with description

    G4UBox(const G4String& pName, G4double pX, G4double pY, G4double pZ);
      // Constructs a box with name, and half lengths pX,pY,pZ

   ~G4UBox();

    void ComputeDimensions(      G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4VSolid* Clone() const;

    inline UBox* GetShape() const;

    inline G4double GetXHalfLength() const;
    inline G4double GetYHalfLength() const;
    inline G4double GetZHalfLength() const;

    inline void SetXHalfLength(G4double dx);
    inline void SetYHalfLength(G4double dy);
    inline void SetZHalfLength(G4double dz);

  public:  // without description

    G4UBox(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UBox(const G4UBox& rhs);
    G4UBox& operator=(const G4UBox& rhs); 
      // Copy constructor and assignment operator.
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline UBox* G4UBox::GetShape() const
{
  return (UBox*) fShape;
}

inline G4double G4UBox::GetXHalfLength() const
{
  return GetShape()->GetXHalfLength();
}
inline G4double G4UBox::GetYHalfLength() const
{
  return GetShape()->GetYHalfLength();
}
inline G4double G4UBox::GetZHalfLength() const
{
  return GetShape()->GetZHalfLength();
}

inline void G4UBox::SetXHalfLength(G4double dx)
{
  GetShape()->SetXHalfLength(dx);
  fRebuildPolyhedron = true;
}
inline void G4UBox::SetYHalfLength(G4double dy)
{
  GetShape()->SetYHalfLength(dy);
  fRebuildPolyhedron = true;
}
inline void G4UBox::SetZHalfLength(G4double dz)
{
  GetShape()->SetZHalfLength(dz);
  fRebuildPolyhedron = true;
}

#endif
