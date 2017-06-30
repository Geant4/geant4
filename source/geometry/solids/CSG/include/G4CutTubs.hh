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
// $Id: G4CutTubs.hh 104316 2017-05-24 13:04:23Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// 
// G4CutTubs is a tube with possible cuts in +-Z.
//           Implementation adapted from G4Tubs (subclass of G4Tubs) and 
//           from TGEo Ctube implementation (by A.Gheata, CERN)
//
// G4CutTubs(pName,pRMin,pRMax,pDZ,pSPhi,pEPhi,pLowNorm,pHighNorm)
//           pName,pRMin,pRMax,pDZ,pSPhi,pEPhi are the same as for G4Tubs,
//           pLowNorm=Outside Normal at -Z
//           pHighNorm=Outsie Normal at +Z.

// Author:   Tatiana Nikitina, CERN
// --------------------------------------------------------------------

#ifndef G4CUTTUBS_HH
#define G4CUTTUBS_HH

#include "G4OTubs.hh"

class G4CutTubs : public G4OTubs
{
  public:  // with description

    G4CutTubs( const G4String& pName,
                     G4double pRMin,
                     G4double pRMax,
                     G4double pDz,
                     G4double pSPhi,
                     G4double pDPhi,
                     G4ThreeVector pLowNorm,
                     G4ThreeVector pHighNorm );
      //
      // Constructs a tubs with the given name and dimensions

   ~G4CutTubs();
      //
      // Destructor

    // Accessors
    
    inline G4ThreeVector GetLowNorm  () const;
    inline G4ThreeVector GetHighNorm () const;  
    
    // Methods for solid

    inline G4double GetCubicVolume();
    inline G4double GetSurfaceArea();

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent( const EAxis pAxis,
                            const G4VoxelLimits& pVoxelLimit,
                            const G4AffineTransform& pTransform,
                                  G4double& pmin, G4double& pmax ) const;

    EInside Inside( const G4ThreeVector& p ) const;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const;

    G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const;
    G4double DistanceToIn(const G4ThreeVector& p) const;
    G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                           const G4bool calcNorm=G4bool(false),
                                 G4bool *validNorm=0, G4ThreeVector *n=0) const;
    G4double DistanceToOut(const G4ThreeVector& p) const;

    G4GeometryType GetEntityType() const;

    G4ThreeVector GetPointOnSurface() const;

    G4VSolid* Clone() const;

    std::ostream& StreamInfo( std::ostream& os ) const;

    // Visualisation functions

    void                DescribeYourselfTo ( G4VGraphicsScene& scene ) const;
    G4Polyhedron*       CreatePolyhedron   () const;

  public:  // without description

    G4CutTubs(__void__&);
      //
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4CutTubs(const G4CutTubs& rhs);
    G4CutTubs& operator=(const G4CutTubs& rhs); 
      // Copy constructor and assignment operator.

  protected:

    G4ThreeVector ApproxSurfaceNormal( const G4ThreeVector& p ) const;
      //
      // Algorithm for SurfaceNormal() following the original
      // specification for points not on the surface

    G4bool IsCrossingCutPlanes() const;
      // Check if the cutted planes are crossing.
      // If 'true' , solid is ill defined
   
    G4double GetCutZ(const G4ThreeVector& p) const;
      // Get Z value of the point on Cutted Plane

    void GetMaxMinZ(G4double& zmin,G4double& zmax)const;
      // Get Max and Min values of Z on Cutted Plane,
      // Used for Calculate BoundingLimits()

  private:

    G4ThreeVector fLowNorm, fHighNorm;
      //
      // Normals of Cut at -/+ Dz

    G4bool fPhiFullCutTube;
      //
      // Flag for identification of section or full tube

    G4double halfCarTolerance, halfRadTolerance, halfAngTolerance;
      //
      // Cached half tolerance values
};

#include "G4CutTubs.icc"

#endif
