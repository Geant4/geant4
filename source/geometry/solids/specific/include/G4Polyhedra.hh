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
// $Id: G4Polyhedra.hh,v 1.22 2010-10-20 08:54:18 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4Polyhedra.hh
//
// Class description:
//
//   Class implementing a CSG-like type "PGON" Geant 3.21 volume,
//   inherited from class G4VCSGfaceted:
//
//   G4Polyhedra( const G4String& name, 
//                G4double phiStart,         - initial phi starting angle
//                G4double phiTotal,         - total phi angle
//                G4int numSide,             - number sides
//                G4int numZPlanes,          - number of z planes
//                const G4double zPlane[],   - position of z planes
//                const G4double rInner[],   - tangent distance to inner surface
//                const G4double rOuter[]  ) - tangent distance to outer surface
//
//   G4Polyhedra( const G4String& name, 
//                G4double phiStart,    - initial phi starting angle
//                G4double phiTotal,    - total phi angle
//                G4int    numSide,     - number sides
//                G4int    numRZ,       - number corners in r,z space
//                const G4double r[],   - r coordinate of these corners
//                const G4double z[] )  - z coordinate of these corners

// Author: 
//   David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#ifndef G4Polyhedra_hh
#define G4Polyhedra_hh

#include "G4VCSGfaceted.hh"
#include "G4PolyhedraSide.hh"

class G4EnclosingCylinder;
class G4ReduciblePolygon;
class G4PolyhedraHistorical
{
  public:

    G4PolyhedraHistorical();
    ~G4PolyhedraHistorical();
    G4PolyhedraHistorical( const G4PolyhedraHistorical &source );
    G4PolyhedraHistorical& operator=( const G4PolyhedraHistorical& right );

    G4double Start_angle;
    G4double Opening_angle;
    G4int   numSide;
    G4int   Num_z_planes;
    G4double *Z_values;
    G4double *Rmin;
    G4double *Rmax;
};

class G4Polyhedra : public G4VCSGfaceted
{
 public:  // with description

  G4Polyhedra( const G4String& name, 
                     G4double phiStart,    // initial phi starting angle
                     G4double phiTotal,    // total phi angle
                     G4int numSide,        // number sides
                     G4int numZPlanes,     // number of z planes
               const G4double zPlane[],    // position of z planes
               const G4double rInner[],    // tangent distance to inner surface
               const G4double rOuter[]  ); // tangent distance to outer surface

  G4Polyhedra( const G4String& name, 
                     G4double phiStart,    // initial phi starting angle
                     G4double phiTotal,    // total phi angle
                     G4int    numSide,     // number sides
                     G4int    numRZ,       // number corners in r,z space
               const G4double r[],         // r coordinate of these corners
               const G4double z[]       ); // z coordinate of these corners

  virtual ~G4Polyhedra();

  // Methods for solid

  EInside Inside( const G4ThreeVector &p ) const;
  G4double DistanceToIn( const G4ThreeVector &p,
                         const G4ThreeVector &v ) const;
  G4double DistanceToIn( const G4ThreeVector &p ) const;
  
  void ComputeDimensions(       G4VPVParameterisation* p,
                          const G4int n,
                          const G4VPhysicalVolume* pRep);

  G4GeometryType  GetEntityType() const;

  G4VSolid* Clone() const;

  G4ThreeVector GetPointOnSurface() const;

  std::ostream& StreamInfo( std::ostream& os ) const;

  G4Polyhedron* CreatePolyhedron() const;
  G4NURBS*      CreateNURBS() const;

  G4bool Reset();

  // Accessors

  inline G4int GetNumSide()     const;
  inline G4double GetStartPhi() const;
  inline G4double GetEndPhi()   const;
  inline G4bool IsOpen()        const;
  inline G4bool IsGeneric()     const;
  inline G4int GetNumRZCorner() const;
  inline G4PolyhedraSideRZ GetCorner( const G4int index ) const;

  inline G4PolyhedraHistorical* GetOriginalParameters() const;
    // Returns internal scaled parameters.
  inline void SetOriginalParameters(G4PolyhedraHistorical* pars);
    // Sets internal parameters. Parameters 'Rmin' and 'Rmax' in input must
    // be scaled first by a factor computed as 'cos(0.5*phiTotal/theNumSide)',
    // if not already scaled.

 public:  // without description

  G4Polyhedra(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

  G4Polyhedra( const G4Polyhedra &source );
  const G4Polyhedra &operator=( const G4Polyhedra &source );
    // Copy constructor and assignment operator.

 protected:  // without description

  inline void SetOriginalParameters();
    // Sets internal parameters for the generic constructor.

  void Create( G4double phiStart,           // initial phi starting angle
               G4double phiTotal,           // total phi angle
               G4int    numSide,            // number sides
               G4ReduciblePolygon *rz );    // rz coordinates
    // Generates the shape and is called by each constructor, after the
    // conversion of the arguments

  void CopyStuff( const G4Polyhedra &source );
  void DeleteStuff();

  // Methods for generation of random points on surface

  G4ThreeVector GetPointOnPlane(G4ThreeVector p0, G4ThreeVector p1,
                                G4ThreeVector p2, G4ThreeVector p3) const;
  G4ThreeVector GetPointOnTriangle(G4ThreeVector p0, G4ThreeVector p1,
                                   G4ThreeVector p2) const;
  G4ThreeVector GetPointOnSurfaceCorners() const;

 protected:  // without description

  G4int   numSide;      // Number of sides
  G4double startPhi;    // Starting phi value (0 < phiStart < 2pi)
  G4double endPhi;      // end phi value (0 < endPhi-phiStart < 2pi)
  G4bool   phiIsOpen;   // true if there is a phi segment
  G4bool   genericPgon; // true if created through the 2nd generic constructor
  G4int   numCorner;    // number RZ points
  G4PolyhedraSideRZ *corners;  // our corners
  G4PolyhedraHistorical  *original_parameters;  // original input parameters

  G4EnclosingCylinder *enclosingCylinder;

};

#include "G4Polyhedra.icc"

#endif
