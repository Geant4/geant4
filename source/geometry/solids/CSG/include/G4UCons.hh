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
// G4UCons
//
// Class description:
//
//   Wrapper class for UCons to make use of UCons from USolids module.

// History:
// 30.10.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------

#ifndef G4UCONS_HH
#define G4UCONS_HH

#include "G4USolid.hh"
#include "UCons.hh"
#include "G4Polyhedron.hh"

class G4UCons : public G4USolid
{
  public:  // with description

    G4UCons(const G4String& pName,
                  G4double pRmin1, G4double pRmax1,
                  G4double pRmin2, G4double pRmax2,
                  G4double pDz,
                  G4double pSPhi, G4double pDPhi);
      // Constructs a cone with the given name and dimensions

   ~G4UCons();

    void ComputeDimensions(       G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep );

    G4VSolid* Clone() const;

    inline UCons* GetShape() const;

    inline G4double GetInnerRadiusMinusZ() const;
    inline G4double GetOuterRadiusMinusZ() const;
    inline G4double GetInnerRadiusPlusZ()  const;
    inline G4double GetOuterRadiusPlusZ()  const;
    inline G4double GetZHalfLength()       const;
    inline G4double GetStartPhiAngle()     const;
    inline G4double GetDeltaPhiAngle()     const;
  
    inline void SetInnerRadiusMinusZ (G4double Rmin1 );
    inline void SetOuterRadiusMinusZ (G4double Rmax1 );
    inline void SetInnerRadiusPlusZ  (G4double Rmin2 );
    inline void SetOuterRadiusPlusZ  (G4double Rmax2 );
    inline void SetZHalfLength       (G4double newDz );
    inline void SetStartPhiAngle     (G4double newSPhi, G4bool trig=true);
    inline void SetDeltaPhiAngle     (G4double newDPhi);

  public:  // without description
       
    G4UCons(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UCons(const G4UCons& rhs);
    G4UCons& operator=(const G4UCons& rhs); 
      // Copy constructor and assignment operator.
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline UCons* G4UCons::GetShape() const
{
  return (UCons*) fShape;
}

inline G4double G4UCons::GetInnerRadiusMinusZ() const
{
  return GetShape()->GetInnerRadiusMinusZ();
}
inline G4double G4UCons::GetOuterRadiusMinusZ() const
{
  return GetShape()->GetOuterRadiusMinusZ();
}
inline G4double G4UCons::GetInnerRadiusPlusZ() const
{
  return GetShape()->GetInnerRadiusPlusZ();
}
inline G4double G4UCons::GetOuterRadiusPlusZ() const
{
  return GetShape()->GetOuterRadiusPlusZ();
}
inline G4double G4UCons::GetZHalfLength() const
{
  return GetShape()->GetZHalfLength();
}
inline G4double G4UCons::GetStartPhiAngle() const
{
  return GetShape()->GetStartPhiAngle();
}
inline G4double G4UCons::GetDeltaPhiAngle() const
{
  return GetShape()->GetDeltaPhiAngle();
}
  
inline void G4UCons::SetInnerRadiusMinusZ(G4double Rmin1)
{
  GetShape()->SetInnerRadiusMinusZ(Rmin1);
  fRebuildPolyhedron = true;
}
inline void G4UCons::SetOuterRadiusMinusZ(G4double Rmax1)
{
  GetShape()->SetOuterRadiusMinusZ(Rmax1);
  fRebuildPolyhedron = true;
}
inline void G4UCons::SetInnerRadiusPlusZ(G4double Rmin2)
{
  GetShape()->SetInnerRadiusPlusZ(Rmin2);
  fRebuildPolyhedron = true;
}
inline void G4UCons::SetOuterRadiusPlusZ(G4double Rmax2)
{
  GetShape()->SetOuterRadiusPlusZ(Rmax2);
  fRebuildPolyhedron = true;
}
inline void G4UCons::SetZHalfLength(G4double newDz)
{
  GetShape()->SetZHalfLength(newDz);
  fRebuildPolyhedron = true;
}
inline void G4UCons::SetStartPhiAngle(G4double newSPhi, G4bool trig)
{
  GetShape()->SetStartPhiAngle(newSPhi, trig);
  fRebuildPolyhedron = true;
}
inline void G4UCons::SetDeltaPhiAngle(G4double newDPhi)
{
  GetShape()->SetDeltaPhiAngle(newDPhi);
  fRebuildPolyhedron = true;
}

#endif
