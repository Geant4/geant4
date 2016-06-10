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
// $Id: G4Sphere.hh 72929 2013-08-14 08:27:27Z gcosmo $
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4USphere
//
// Class description:
//
//   Wrapper class for USphere to make use of USphere from USolids module.

// History:
// 13.09.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------
#ifndef G4USPHERE_HH
#define G4USPHERE_HH

#include "G4USolid.hh"
#include "USphere.hh"
#include "G4Polyhedron.hh"

class G4USphere : public G4USolid
{
  public:  // with description

    G4USphere(const G4String& pName,
                    G4double pRmin, G4double pRmax,
                    G4double pSPhi, G4double pDPhi,
                    G4double pSTheta, G4double pDTheta);
      // Constructs a sphere or sphere shell section
      // with the given name and dimensions
       
   ~G4USphere();

    void ComputeDimensions(      G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4VSolid* Clone() const;

    inline USphere* GetShape() const;
       
    inline G4double GetInnerRadius    () const;
    inline G4double GetOuterRadius    () const;
    inline G4double GetStartPhiAngle  () const;
    inline G4double GetDeltaPhiAngle  () const;
    inline G4double GetStartThetaAngle() const;
    inline G4double GetDeltaThetaAngle() const;

    inline void SetInnerRadius    (G4double newRMin);
    inline void SetOuterRadius    (G4double newRmax);
    inline void SetStartPhiAngle  (G4double newSphi, G4bool trig=true);
    inline void SetDeltaPhiAngle  (G4double newDphi);
    inline void SetStartThetaAngle(G4double newSTheta);
    inline void SetDeltaThetaAngle(G4double newDTheta);

  public:  // without description
   
    G4USphere(__void__&);
      //
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4USphere(const G4USphere& rhs);
    G4USphere& operator=(const G4USphere& rhs); 
      // Copy constructor and assignment operator.
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline USphere* G4USphere::GetShape() const
{
  return (USphere*) fShape;
}

inline G4double G4USphere::GetInnerRadius() const
{
  return GetShape()->GetInnerRadius();
}
inline G4double G4USphere::GetOuterRadius() const
{
  return GetShape()->GetOuterRadius();
}
inline G4double G4USphere::GetStartPhiAngle() const
{
  return GetShape()->GetStartPhiAngle();
}
inline G4double G4USphere::GetDeltaPhiAngle() const
{
  return GetShape()->GetDeltaPhiAngle();
}
inline G4double G4USphere::GetStartThetaAngle() const
{
  return GetShape()->GetStartThetaAngle();
}
inline G4double G4USphere::GetDeltaThetaAngle() const
{
  return GetShape()->GetDeltaThetaAngle();
}

inline void G4USphere::SetInnerRadius(G4double newRMin)
{
  GetShape()->SetInnerRadius(newRMin);
  fRebuildPolyhedron = true;
}
inline void G4USphere::SetOuterRadius(G4double newRmax)
{
  GetShape()->SetOuterRadius(newRmax);
  fRebuildPolyhedron = true;
}
inline void G4USphere::SetStartPhiAngle(G4double newSphi, G4bool trig)
{
  GetShape()->SetStartPhiAngle(newSphi, trig);
  fRebuildPolyhedron = true;
}
inline void G4USphere::SetDeltaPhiAngle(G4double newDphi)
{
  GetShape()->SetDeltaPhiAngle(newDphi);
  fRebuildPolyhedron = true;
}
inline void G4USphere::SetStartThetaAngle(G4double newSTheta)
{
  GetShape()->SetStartThetaAngle(newSTheta);
  fRebuildPolyhedron = true;
}
inline void G4USphere::SetDeltaThetaAngle(G4double newDTheta)
{
  GetShape()->SetDeltaThetaAngle(newDTheta);
  fRebuildPolyhedron = true;
}

#endif
