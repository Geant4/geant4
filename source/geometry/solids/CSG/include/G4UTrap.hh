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
// G4UTrap
//
// Class description:
//
//   Wrapper class for UTrap to make use of UTrap from USolids module.

// History:
// 13.09.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------
#ifndef G4UTrap_HH
#define G4UTrap_HH

#include "G4USolid.hh"
#include "UTrap.hh"

class G4Polyhedron;

class G4UTrap : public G4USolid 
{
  public:  // with description

    G4UTrap( const G4String& pName,
                   G4double pDz,
                   G4double pTheta, G4double pPhi,
                   G4double pDy1, G4double pDx1, G4double pDx2,
                   G4double pAlp1,
                   G4double pDy2, G4double pDx3, G4double pDx4,
                   G4double pAlp2 );
      //
      // The most general constructor for G4Trap which prepares plane
      // equations and corner coordinates from parameters

    G4UTrap( const G4String& pName,
             const G4ThreeVector pt[8] ) ;
      //
      // Prepares plane equations and parameters from corner coordinates

    G4UTrap( const G4String& pName,
                   G4double pZ,
                   G4double pY,
                   G4double pX, G4double pLTX );
      //
      // Constructor for Right Angular Wedge from STEP (assumes pLTX<=pX)

    G4UTrap( const G4String& pName,
                   G4double pDx1,  G4double pDx2,
                   G4double pDy1,  G4double pDy2,
                   G4double pDz );
      //
      // Constructor for G4Trd       

    G4UTrap(const G4String& pName,
                  G4double pDx, G4double pDy, G4double pDz,
                  G4double pAlpha, G4double pTheta, G4double pPhi );
      //
      // Constructor for G4Para

    G4UTrap( const G4String& pName );
      //
      // Constructor for "nominal" G4Trap whose parameters are to be set
      // by a G4VPVParamaterisation later

   ~G4UTrap();

    void ComputeDimensions(      G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4VSolid* Clone() const;

    inline UTrap* GetShape() const;

    inline G4double GetZHalfLength()  const;
    inline G4double GetYHalfLength1() const;
    inline G4double GetXHalfLength1() const;
    inline G4double GetXHalfLength2() const;
    inline G4double GetTanAlpha1()    const;
    inline G4double GetYHalfLength2() const;
    inline G4double GetXHalfLength3() const;
    inline G4double GetXHalfLength4() const;
    inline G4double GetTanAlpha2()    const;
    inline TrapSidePlane GetSidePlane(G4int n) const;
    inline G4ThreeVector GetSymAxis() const;

    inline void SetAllParameters(G4double pDz, G4double pTheta, G4double pPhi,
                                 G4double pDy1, G4double pDx1, G4double pDx2,
                                 G4double pAlp1,
                                 G4double pDy2, G4double pDx3, G4double pDx4,
                                 G4double pAlp2);
    inline void SetPlanes(const G4ThreeVector pt[8]);

  public:  // without description

    G4UTrap(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UTrap(const G4UTrap& rhs);
    G4UTrap& operator=(const G4UTrap& rhs); 
      // Copy constructor and assignment operator.

    G4Polyhedron* CreatePolyhedron   () const;
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline UTrap* G4UTrap::GetShape() const
{
  return (UTrap*) fShape;
}

inline G4double G4UTrap::GetZHalfLength() const
{
  return GetShape()->GetZHalfLength();
}
inline G4double G4UTrap::GetYHalfLength1() const
{
  return GetShape()->GetYHalfLength1();
}
inline G4double G4UTrap::GetXHalfLength1() const
{
  return GetShape()->GetXHalfLength1();
}
inline G4double G4UTrap::GetXHalfLength2() const
{
  return GetShape()->GetXHalfLength2();
}
inline G4double G4UTrap::GetTanAlpha1() const
{
  return GetShape()->GetTanAlpha1();
}
inline G4double G4UTrap::GetYHalfLength2() const
{
  return GetShape()->GetYHalfLength2();
}
inline G4double G4UTrap::GetXHalfLength3() const
{
  return GetShape()->GetXHalfLength3();
}
inline G4double G4UTrap::GetXHalfLength4() const
{
  return GetShape()->GetXHalfLength4();
}
inline G4double G4UTrap::GetTanAlpha2() const
{
  return GetShape()->GetTanAlpha2();
}
inline TrapSidePlane G4UTrap::GetSidePlane(G4int n) const
{
  UTrapSidePlane iplane = GetShape()->GetSidePlane(n);
  TrapSidePlane oplane = {iplane.a, iplane.b, iplane.c, iplane.d };
  return oplane;
}
inline G4ThreeVector G4UTrap::GetSymAxis() const
{
  UVector3 axis = GetShape()->GetSymAxis();
  return G4ThreeVector(axis.x(), axis.y(), axis.z());
}

inline
void G4UTrap::SetAllParameters(G4double pDz, G4double pTheta, G4double pPhi,
                               G4double pDy1, G4double pDx1, G4double pDx2,
                               G4double pAlp1,
                               G4double pDy2, G4double pDx3, G4double pDx4,
                               G4double pAlp2)
{
  GetShape()->SetAllParameters(pDz, pTheta, pPhi,
                               pDy1, pDx1, pDx2, pAlp1,
                               pDy2, pDx3, pDx4, pAlp2);
  fRebuildPolyhedron = true;
}

inline void G4UTrap::SetPlanes(const G4ThreeVector pt[8])
{
  UVector3 upt[8];
  for (unsigned int i=0; i<8; ++i)
  {
    upt[i] = UVector3(pt[i].x(), pt[i].y(), pt[i].z());
  }
  GetShape()->SetPlanes(upt);
  fRebuildPolyhedron = true;
}

#endif
