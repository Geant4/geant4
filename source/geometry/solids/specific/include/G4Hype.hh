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
// $Id: G4Hype.hh 104316 2017-05-24 13:04:23Z gcosmo $
// $Original: G4Hype.hh,v 1.0 1998/06/09 16:57:50 safai Exp $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4Hype
//
// Class description:
//
//   This class implements a tube with hyperbolic profile.
//
//   It describes an hyperbolic volume with curved sides parallel to
//   the z-axis. The solid has a specified half-length along the z axis,
//   about which it is centered, and a given minimum and maximum radius.
//   A minimum radius of 0 signifies a filled Hype (with hyperbolical
//   inner surface). To have a filled Hype the user must specify 
//   inner radius = 0 AND inner stereo angle = 0.
// 
//   The inner and outer hyperbolical surfaces can have different
//   stereo angles. A stereo angle of 0 gives a cylindrical surface.

// Authors: 
//      Ernesto Lamanna (Ernesto.Lamanna@roma1.infn.it) &
//      Francesco Safai Tehrani (Francesco.SafaiTehrani@roma1.infn.it)
//      Rome, INFN & University of Rome "La Sapienza",  9 June 1998.
//
// --------------------------------------------------------------------
#ifndef G4HYPE_HH
#define G4HYPE_HH

#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Polyhedron.hh"

class G4SolidExtentList;
class G4ClippablePolygon;

class G4Hype : public G4VSolid
{
 public:  // with description

  G4Hype(const G4String& pName,
               G4double  newInnerRadius,
               G4double  newOuterRadius,
               G4double  newInnerStereo,
               G4double  newOuterStereo,
               G4double  newHalfLenZ);

  virtual ~G4Hype();
    
  void ComputeDimensions(G4VPVParameterisation* p,
                         const G4int n,
                         const G4VPhysicalVolume* pRep);

  void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

  G4bool CalculateExtent(const EAxis pAxis,
                         const G4VoxelLimits& pVoxelLimit,
                         const G4AffineTransform& pTransform,
                               G4double& pMin, G4double& pMax) const;

  inline G4double GetInnerRadius () const;
  inline G4double GetOuterRadius () const;
  inline G4double GetZHalfLength () const;
  inline G4double GetInnerStereo () const;
  inline G4double GetOuterStereo () const;

  inline void SetInnerRadius (G4double newIRad);
  inline void SetOuterRadius (G4double newORad);
  inline void SetZHalfLength (G4double newHLZ);
  inline void SetInnerStereo (G4double newISte);
  inline void SetOuterStereo (G4double newOSte);

  EInside Inside(const G4ThreeVector& p) const;

  G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const;

  G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const;
  G4double DistanceToIn(const G4ThreeVector& p) const;
  G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                         const G4bool calcNorm=G4bool(false),
                         G4bool *validNorm=0, G4ThreeVector *n=0) const;
  G4double DistanceToOut(const G4ThreeVector& p) const;

  G4GeometryType  GetEntityType() const;

  G4VSolid* Clone() const;

  std::ostream& StreamInfo(std::ostream& os) const;
  
  G4double GetCubicVolume();
  G4double GetSurfaceArea();

  G4ThreeVector GetPointOnSurface() const;

  void          DescribeYourselfTo (G4VGraphicsScene& scene) const;
  G4VisExtent   GetExtent          () const;
  G4Polyhedron* CreatePolyhedron   () const;
  G4Polyhedron* GetPolyhedron      () const;

 public:  // without description

  G4Hype(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

  G4Hype(const G4Hype& rhs);
  G4Hype& operator=(const G4Hype& rhs); 
    // Copy constructor and assignment operator.

 protected:  // without description
  
  inline G4bool InnerSurfaceExists() const;
    // whether we have an inner surface or not

  static G4double ApproxDistOutside( G4double pr, G4double pz,
                                     G4double r0, G4double tanPhi );
  static G4double ApproxDistInside( G4double pr, G4double pz,
                                    G4double r0, G4double tan2Phi );
    // approximate isotropic distance to hyperbolic surface 

  inline G4double HypeInnerRadius2(G4double zVal) const;
  inline G4double HypeOuterRadius2(G4double zVal) const;
    // values of hype radius at a given Z

  static G4int IntersectHype( const G4ThreeVector &p, const G4ThreeVector &v, 
                              G4double r2, G4double tan2Phi, G4double s[2] );
    // intersection with hyperbolic surface

 protected:

  G4double innerRadius;
  G4double outerRadius;
  G4double halfLenZ;
  G4double innerStereo;
  G4double outerStereo;

  // precalculated parameters, squared quantities

  G4double tanInnerStereo;
  G4double tanOuterStereo;
  G4double tanInnerStereo2; // squared tan of Inner Stereo angle
  G4double tanOuterStereo2; // squared tan of Outer Stereo angle
  G4double innerRadius2;    // squared Inner Radius
  G4double outerRadius2;    // squared Outer Radius
  G4double endInnerRadius2; // squared endcap Inner Radius
  G4double endOuterRadius2; // squared endcap Outer Radius
  G4double endInnerRadius; // endcap Inner Radius
  G4double endOuterRadius; // endcap Outer Radius
  
  // Used by distanceToOut

  enum ESide {outerFace,innerFace,leftCap, rightCap};

 private:

  G4double asinh(G4double arg);

 private:
  
  G4double fCubicVolume;
  G4double fSurfaceArea;

  G4double fHalfTol;

  mutable G4bool fRebuildPolyhedron;
  mutable G4Polyhedron* fpPolyhedron;
};

#include "G4Hype.icc"

#endif
