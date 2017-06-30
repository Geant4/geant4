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
// $Id: G4GenerciPolycone.hh 72077 2013-07-08 12:31:44Z tnikitin $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4GenericPolycone
//
// Class description:
//
//   Class implementing a GenericPolycone constructed by points with
//          (r,z)coordinates, allows Z 'go back' 
//
//
//   G4GenericPolycone( const G4String& name, 
//               G4double phiStart,   // initial phi starting angle
//               G4double phiTotal,   // total phi angle
//               G4int    numRZ,      // number corners in r,z space
//               const G4double r[],  // r coordinate of these corners
//               const G4double z[])  // z coordinate of these corners
//
//  
// --------------------------------------------------------------------
#ifndef G4GenericPolycone_hh
#define G4GenericPolycone_hh

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UGENERICPOLYCONE 1
#endif

#if defined(G4GEOM_USE_UGENERICPOLYCONE)
  #define G4UGenericPolycone G4GenericPolycone
  #include "G4UGenericPolycone.hh"
#else

#include "G4VCSGfaceted.hh"
#include "G4PolyconeSide.hh"

class G4EnclosingCylinder;
class G4ReduciblePolygon;
class G4VCSGface;

class G4GenericPolycone : public G4VCSGfaceted 
{

 public:  // with description

  G4GenericPolycone( const G4String& name, 
                    G4double phiStart,    // initial phi starting angle
                    G4double phiTotal,    // total phi angle
                    G4int    numRZ,       // number corners in r,z space
                    const G4double r[],         // r coordinate of these corners
                    const G4double z[]       ); // z coordinate of these corners

  virtual ~G4GenericPolycone();
  
  // Methods for solid

  EInside Inside( const G4ThreeVector &p ) const;
  G4double DistanceToIn( const G4ThreeVector &p, const G4ThreeVector &v ) const;
  G4double DistanceToIn( const G4ThreeVector &p ) const;

  void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;
  G4bool CalculateExtent(const EAxis pAxis,
                         const G4VoxelLimits& pVoxelLimit,
                         const G4AffineTransform& pTransform,
                         G4double& pmin, G4double& pmax) const;

  G4ThreeVector GetPointOnSurface() const;

  G4GeometryType GetEntityType() const;

  G4VSolid* Clone() const;

  std::ostream& StreamInfo(std::ostream& os) const;

  G4Polyhedron* CreatePolyhedron() const;

  G4bool Reset();

  // Accessors

  inline G4double GetStartPhi()    const;
  inline G4double GetEndPhi()      const;
  inline G4double GetSinStartPhi() const;
  inline G4double GetCosStartPhi() const;
  inline G4double GetSinEndPhi()   const;
  inline G4double GetCosEndPhi()   const;
  inline G4bool IsOpen()           const;
  inline G4int  GetNumRZCorner()   const;
  inline G4PolyconeSideRZ GetCorner(G4int index) const;
  
 public:  // without description

  G4GenericPolycone(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

  G4GenericPolycone( const G4GenericPolycone &source );
  G4GenericPolycone &operator=( const G4GenericPolycone &source );
    // Copy constructor and assignment operator.

 protected:  // without description

  // Generic initializer, called by all constructors

  // void SetOriginalParameters(G4ReduciblePolygon *rz);

  void Create( G4double phiStart,        // initial phi starting angle
               G4double phiTotal,        // total phi angle
               G4ReduciblePolygon *rz ); // r/z coordinate of these corners

  void CopyStuff( const G4GenericPolycone &source );

  // Methods for random point generation

 protected:  // without description

  // Here are our parameters

  G4double startPhi;    // Starting phi value (0 < phiStart < 2pi)
  G4double endPhi;      // end phi value (0 < endPhi-phiStart < 2pi)
  G4bool   phiIsOpen;   // true if there is a phi segment
  G4int    numCorner;   // number RZ points
  G4PolyconeSideRZ *corners;  // corner r,z points
 
  // Our quick test

  G4EnclosingCylinder *enclosingCylinder;
  
};

#include "G4GenericPolycone.icc"

#endif

#endif
