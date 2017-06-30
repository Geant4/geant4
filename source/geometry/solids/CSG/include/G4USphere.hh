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

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

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

    G4double GetInnerRadius    () const;
    G4double GetOuterRadius    () const;
    G4double GetStartPhiAngle  () const;
    G4double GetDeltaPhiAngle  () const;
    G4double GetStartThetaAngle() const;
    G4double GetDeltaThetaAngle() const;
    G4double GetSinStartPhi    () const;
    G4double GetCosStartPhi    () const;
    G4double GetSinEndPhi      () const;
    G4double GetCosEndPhi      () const;
    G4double GetSinStartTheta  () const;
    G4double GetCosStartTheta  () const;
    G4double GetSinEndTheta    () const;
    G4double GetCosEndTheta    () const;

    void SetInnerRadius    (G4double newRMin);
    void SetOuterRadius    (G4double newRmax);
    void SetStartPhiAngle  (G4double newSphi, G4bool trig=true);
    void SetDeltaPhiAngle  (G4double newDphi);
    void SetStartThetaAngle(G4double newSTheta);
    void SetDeltaThetaAngle(G4double newDTheta);

    inline G4GeometryType GetEntityType() const;

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const;

    G4Polyhedron* CreatePolyhedron() const;

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

inline G4GeometryType G4USphere::GetEntityType() const
{
  return "G4Sphere";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
